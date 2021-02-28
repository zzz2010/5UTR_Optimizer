#!/usr/bin/env Rscript

#Rscript evaldModel_final.R -p output/input.fa -a data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt -k 5 -r 0.1 -m 1,2 -o output/muscle_randomforest.model
library("optparse")
library(e1071)
library(glmnet)
library(class)
library(randomForest)
library(rpart)
library(caret)
library(methods)
library(Matrix)
library(foreach)
library(doMC)
library(gridExtra)
library(Metrics)
library(xlsx)
 
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="sequence feature files prefix", metavar="character"),
    make_option(c("-k", "--min_rna_rpkm"), default=50.0, type="double",
              help="minimal RNASEQ RPKM for the gene included in training. In paper, we used 5 for muscle cell", metavar="min_rna_rpkm"),

    make_option(c("-r", "--min_riboseq_rpkm"), default=5.0, type="double",
              help="minimal RiboSEQ RPKM for the gene included in training. In paper, we used 0.1 for muscle cell", metavar="min_riboseq_rpkm"),

    make_option(c("-a", "--annotation_file"), type="character", default=NULL, 
              help="transcript annotation file with gene id and RPKM information for rnaseq and riboseq", metavar="character"),

    make_option(c("-m", "--modellist"), type="character", default="1,2,3,4",
              help="1: randomforest, 2: glmnet, 3: regression tree, 4: SVM", metavar="modellist"),

    make_option(c("-j", "--n_cpus"), type="integer", default=16,
              help="number of cpus", metavar="n_cpus"),
    make_option(c("-o", "--out"), type="character", default="out.xlsx",
              help="output evaluation result excel file", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

registerDoMC(cores=opt$n_cpus)

prefix=opt$prefix

annotFile=opt$annotation_file

modellist=opt$modellist


feat.df=read.table(sprintf("%s.sparseFeature.txt.gz",prefix))
feat.mat=sparseMatrix(i=feat.df[,1],j=feat.df[,2],x=feat.df[,3],index1=F)
rownames(feat.mat)=read.table(sprintf("%s.sparseFeature.rowname",prefix),row.names=1)[,1]
colnames(feat.mat)=read.table(sprintf("%s.sparseFeature.colname",prefix),row.names=1)[1:ncol(feat.mat),1] ##because last few Feature(kmer) may not occurs in the data, so need to trim



TE.df=read.table(annotFile,row.names=1,header=T)

mRNA.RPKM.filter=TE.df[,'rpkm_rnaseq']> opt$min_rna_rpkm
ribo.RPKM.filter=TE.df[,'rpkm_riboseq']> opt$min_riboseq_rpkm

nomiss=complete.cases(TE.df)
TE.df=TE.df[rownames(TE.df) %in% rownames(feat.mat) & mRNA.RPKM.filter  & ribo.RPKM.filter & nomiss,]
if(nrow(TE.df)<30){
  print("Too few data passing the filter criteria, please relax the filter criteria and try again !")
}

sel.row=match(row.names(TE.df),rownames(feat.mat))

selfeat.mat=feat.mat[sel.row,]
te=TE.df[,'te']
ribo=TE.df[,'rpkm_riboseq']
rnaseq=TE.df[,'rpkm_rnaseq']

print(sprintf("Number of training samples: %d",length(te)))


predictive_performance<-function(y){

test.id.list=createFolds(y,k=10)

featRange=1:ncol(selfeat.mat) #31

########Cross validation to evaluate different algorithm performance ##########
performances=foreach(i=1:length(test.id.list))%dopar%{
test.id=unlist(test.id.list[i])
train.x=selfeat.mat[-test.id,featRange]
train.y=y[-test.id]
test.x=selfeat.mat[test.id,featRange]
test.y=y[test.id]

res_list=foreach (model = unlist(strsplit(modellist, ",")))%do%{
    model=as.integer(model)
    if(model==1){
    #########random forest#######
    fit=randomForest(as.matrix(train.x), train.y,importance=T)
    pred.y=predict(fit,test.x)
    model_name="randomforest"
    }


    if(model==2){
      nfolds=5
    #########glmnet##############
    ##train a model
      fit=cv.glmnet(train.x,train.y,nfolds=nfolds) #,parallel=F,alpha=0
      pred.y=predict(fit$glmnet.fit,newx=test.x,s=fit$lambda.min)
      model_name="glmnet"
    }

    if(model==3){
    #######regression tree#######
    df=data.frame(y=train.y,as.matrix(train.x))
    fit=prune(rpart(y ~ .,df),cp=0.01)
    pred.y=predict(fit,newdata=data.frame(y=NA,as.matrix(test.x)))
    model_name="rpart"
    }

    if(model==4){
      fit=svm(train.x,train.y)
      pred.y=predict(fit,test.x)
      model_name="svm"
    }

    ##evaluate
    sp_cor=cor(pred.y,test.y,method="spearman")
    rsq=cor(pred.y,test.y)^2
   print(c(model_name,sp_cor))
   data.frame(sp_cor=sp_cor,rsq=rsq,model=model_name,cv_k=i)

    }
  do.call(rbind, res_list)
  }
do.call(rbind, performances)
}



### build full model and save ###
###model translation efficiency

#plotTable<-function (d,title){
#table <- tableGrob(d)
#
#grid.newpage()
#h <- grobHeight(table)
#w <- grobWidth(table)
#title <- textGrob(title, y=unit(0.5,"npc") + 0.5*h,
#                  vjust=0, gp=gpar(fontsize=20))
#footnote <- textGrob("footnote",
#                     x=unit(0.5,"npc") - 0.5*w,
#                     y=unit(0.5,"npc") - 0.5*h,
#                  vjust=1, hjust=0,gp=gpar( fontface="italic"))
#gt <- gTree(children=gList(table, title, footnote))
#grid.draw(gt)
#}








print("building TE model.....")
y=log(te)
evalres_te=(predictive_performance(y))
write.xlsx(evalres_te, file=opt$out, sheetName="TE_result", row.names=FALSE)
full.model.te=randomForest(as.matrix(selfeat.mat),log(te),importance=T)
#print("Top features for predicting Translation Efficiency")
featScore=importance(full.model.te, type=1)
write.xlsx( as.data.frame(t(featScore)), file=opt$out, sheetName="TE_features_accuracy",append=TRUE, row.names=FALSE)

featScore=importance(full.model.te, type=1)
write.xlsx( as.data.frame(t(featScore)), file=opt$out, sheetName="TE_features_impurity",append=TRUE, row.names=FALSE)
#featScore[order(-featScore)[1:10],]
#featScore=importance(full.model.te, type=2)
#featScore[order(-featScore)[1:10],]

###model gene expression
print("building RNA expression model.....")

y=log(rnaseq)
evalres_rna=(predictive_performance(y))
write.xlsx(evalres_rna, file=opt$out, sheetName="RNA_result", append=TRUE, row.names=FALSE)

full.model.rna=randomForest(as.matrix(selfeat.mat),log(rnaseq),importance=T)

#print("Top features for predicting mRNA expression")
featScore=importance(full.model.rna, type=1)
write.xlsx( as.data.frame(t(featScore)), file=opt$out, sheetName="RNA_features_accuracy",append=TRUE, row.names=FALSE)

featScore=importance(full.model.rna, type=1)
write.xlsx( as.data.frame(t(featScore)), file=opt$out, sheetName="RNA_features_impurity",append=TRUE, row.names=FALSE)
#featScore[order(-featScore)[1:10],]
#featScore=importance(full.model.rna, type=2)
#featScore[order(-featScore)[1:10],]


