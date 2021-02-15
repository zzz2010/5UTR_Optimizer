#!/usr/bin/env Rscript

#Rscript buildModel_final.R -p output/input.fa -a data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt -k 5 -r 0.1 -m 1 -o output/muscle_randomforest.model
library("optparse")
library(methods)
library(Matrix)
library(foreach)
library(doMC)
library(Metrics)
library(e1071)
library(glmnet)
library(class)
library(randomForest)
library(rpart)
library(caret)
 
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="sequence feature files prefix", metavar="character"),
    make_option(c("-k", "--min_rna_rpkm"), default=50,
              help="minimal RNASEQ RPKM for the gene included in training. In paper, we used 5 for muscle cell", metavar="min_rna_rpkm"),

    make_option(c("-r", "--min_riboseq_rpkm"), default=5,
              help="minimal RiboSEQ RPKM for the gene included in training. In paper, we used 0.1 for muscle cell", metavar="min_ribo_rpkm"),

    make_option(c("-a", "--annotation_file"), type="character", default=NULL, 
              help="transcript annotation file with gene id and RPKM information for rnaseq and riboseq", metavar="character"),

    make_option(c("-m", "--model"), type="integer", default=1,
              help="1: randomforest, 2: glmnet, 3: regression tree, 4: SVM", metavar="model"),
    make_option(c("-j", "--n_cpus"), type="integer", default=16,
              help="number of cpus", metavar="n_cpus"),
    make_option(c("-o", "--out"), type="character", default="out.model",
              help="output file name [default= %default]", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

registerDoMC(cores=opt$n_cpus)

prefix=opt$prefix

annotFile=opt$annotation_file

model=opt$model


feat.df=read.table(sprintf("%s.sparseFeature.txt.gz",prefix))
feat.mat=sparseMatrix(i=feat.df[,1],j=feat.df[,2],x=feat.df[,3],index1=F)
rownames(feat.mat)=read.table(sprintf("%s.sparseFeature.rowname",prefix),row.names=1)[,1]
colnames(feat.mat)=read.table(sprintf("%s.sparseFeature.colname",prefix),row.names=1)[1:ncol(feat.mat),1] ##because last few Feature(kmer) may not occurs in the data, so need to trim

#cell=".HEK_Andrev2015" ##HEK cell
#cell=".pc3"
#cell=""  ## muscle  

TE.df=read.table(annotFile,row.names=1,header=T)
#TE.df=read.table("data/df_counts_and_len.TE_sorted.with_annot.txt",row.names=1,header=T)

#TE.df=read.table(sprintf("data/df_counts_and_len.TE_sorted%s.with_annot.txt",cell),row.names=1,header=T)

mRNA.RPKM.filter=TE.df[,'rpkm_rnaseq']>opt$min_rna_rpkm
ribo.RPKM.filter=TE.df[,'rpkm_riboseq']>opt$min_riboseq_rpkm

nomiss=complete.cases(TE.df)
TE.df=TE.df[rownames(TE.df) %in% rownames(feat.mat) & mRNA.RPKM.filter & nomiss,]

sel.row=match(row.names(TE.df),rownames(feat.mat))

selfeat.mat=feat.mat[sel.row,]
te=TE.df[,'te']
ribo=TE.df[,'rpkm_riboseq']
rnaseq=TE.df[,'rpkm_rnaseq']

print(sprintf("Number of training samples: %d",length(te)))



build_model<-function(y){

featRange=1:ncol(selfeat.mat) #31
train.x=selfeat.mat[,featRange]
train.y=y[-test.id]

if(model==1){
#########glmnet##############
##train a model
fit=cv.glmnet(train.x,train.y,nfolds=nfolds,parallel=T,alpha=0)

}
if(model==2){
#########random forest#######
fit=randomForest(as.matrix(train.x), train.y,importance=T)

}
if(model==3){
#######regression tree#######
df=data.frame(y=train.y,as.matrix(train.x))
fit=prune(rpart(y ~ .,df),cp=0.01)

}
if(model==4){
fit=svm(train.x,train.y)	

}

return(fit)

}


### build full model and save ###
###model translation efficiency
print("building TE model.....")
y=log(te)
full.model.te=build_model(y)


###model gene expression
print("building RNA expression model.....")

y=log(rnaseq)
full.model.rna=build_model(y)


modelFn=opt$out
save(full.model.te,full.model.rna,file=modelFn)
print(sprintf("save model to file : %s",modelFn))
