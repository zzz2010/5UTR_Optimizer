library(methods)
library(Matrix)
library(foreach)
library(doMC)
library(Metrics)
registerDoMC(cores=30)
prefix="output/gencode_v17_5utr_15bpcds.fa"
feat.df=read.table(sprintf("%s.sparseFeature.txt.gz",prefix))
feat.mat=sparseMatrix(i=feat.df[,1],j=feat.df[,2],x=feat.df[,3],index1=F)
rownames(feat.mat)=read.table(sprintf("%s.sparseFeature.rowname",prefix),row.names=1)[,1]
colnames(feat.mat)=read.table(sprintf("%s.sparseFeature.colname",prefix),row.names=1)[1:ncol(feat.mat),1] ##because last few Feature(kmer) may not occurs in the data, so need to trim

cell=".HEK_Andrev2015" ##HEK cell
#cell=".pc3"
#cell=""  ## muscle  

#TE.df=read.table("data/df_counts_and_len.TE_sorted.with_annot.txt",row.names=1,header=T)

TE.df=read.table(sprintf("data/df_counts_and_len.TE_sorted%s.with_annot.txt",cell),row.names=1,header=T)

if(cell==""){
mRNA.RPKM.filter=TE.df[,'rpkm_rnaseq']>5
ribo.RPKM.filter=TE.df[,'rpkm_riboseq']>0.1
}else{
mRNA.RPKM.filter=TE.df[,'rpkm_rnaseq']>50
ribo.RPKM.filter=TE.df[,'rpkm_riboseq']>5
}
nomiss=complete.cases(TE.df)
TE.df=TE.df[rownames(TE.df) %in% rownames(feat.mat) & mRNA.RPKM.filter & nomiss,]

sel.row=match(row.names(TE.df),rownames(feat.mat))

selfeat.mat=feat.mat[sel.row,]
te=TE.df[,'te']
ribo=TE.df[,'rpkm_riboseq']
rnaseq=TE.df[,'rpkm_rnaseq']

print(sprintf("Number of training samples: %d",length(te)))
library(e1071)
library(glmnet)
library(class)
library(randomForest)
library(rpart)
library(caret)


predictive_performance<-function(y){
model=2  ##random forest seems the best
#model=1 ## glmnet
test.id.list=createFolds(y,k=10)

featRange=1:ncol(selfeat.mat) #31

########Cross validation to evaluate different algorithm performance ##########
performances=foreach(i=1:length(test.id.list))%dopar%{
test.id=unlist(test.id.list[i])
train.x=selfeat.mat[-test.id,featRange]
train.y=y[-test.id]
test.x=selfeat.mat[test.id,featRange]
test.y=y[test.id]

nfolds=5
if(model==1){
#########glmnet##############
##train a model
fit=cv.glmnet(train.x,train.y,nfolds=nfolds,parallel=T,alpha=0)
pred.y=predict(fit$glmnet.fit,newx=test.x,s=fit$lambda.min)
}
if(model==2){
#########random forest#######
fit=randomForest(as.matrix(train.x), train.y)
pred.y=predict(fit,test.x)
}
if(model==3){
#######regression tree#######
df=data.frame(y=train.y,as.matrix(train.x))
fit=prune(rpart(y ~ .,df),cp=0.01)
pred.y=predict(fit,newdata=data.frame(y=NA,as.matrix(test.x)))
}
if(model==4){
fit=svm(train.x,train.y)	
pred.y=predict(fit,test.x)
}

##evaluate
#cor(pred.y,test.y,method="spearman")
cor(pred.y,test.y)
}
print(mean(na.omit(unlist(performances))))
print(length(y))
}



### build full model and save ###
###model translation efficiency
print("building TE model.....")
y=log(te)
predictive_performance(y)

full.model.te=randomForest(as.matrix(selfeat.mat),log(te),importance=T) 
print("Top features for predicting Translation Efficiency")
featScore=importance(full.model.te, type=1)
featScore[order(-featScore)[1:10],]
featScore=importance(full.model.te, type=2)
featScore[order(-featScore)[1:10],]

###model gene expression
print("building RNA expression model.....")

y=log(rnaseq)
predictive_performance(y)

full.model.rna=randomForest(as.matrix(selfeat.mat),log(rnaseq),importance=T)
print("Top features for predicting mRNA expression")
featScore=importance(full.model.rna, type=1)
featScore[order(-featScore)[1:10],]
featScore=importance(full.model.rna, type=2)
featScore[order(-featScore)[1:10],]

modelFn=sprintf("%s%s.big.model",prefix,cell)
save(full.model.te,full.model.rna,file=modelFn)
print(sprintf("save model to file : %s",modelFn))
