###Author: Zhizhuo Zhang
###Date: 12/27/2015
###use bigger training data 
###smaller pop size, and fewer generations
###run more rounds

#Rscript evolutionDesign.R -p output/input.fa -a data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt -t ribo -k 5 -r 0.1 -m output/muscle_randomforest.model  -o output/
library(GA)
library("optparse")
library(randomForest)
library(methods)
library(Matrix)
library(foreach)
library(doMC)
library(Metrics)
library(doParallel)
library(parallel)
library(seqinr)
library(Biostrings)
options(stringsAsFactors = FALSE)

option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="sequence feature files prefix", metavar="character"),
    make_option(c("-k", "--min_rna_rpkm"), default=50,
              help="minimal RNASEQ RPKM for the gene included in training. In paper, we used 5 for muscle cell", metavar="min_rna_rpkm"),

    make_option(c("-r", "--min_riboseq_rpkm"), default=5,
              help="minimal RiboSEQ RPKM for the gene included in training. In paper, we used 0.1 for muscle cell", metavar="min_riboseq_rpkm"),

    make_option(c("-a", "--annotation_file"), type="character", default=NULL,
              help="transcript annotation file with gene id and RPKM information for rnaseq and riboseq", metavar="character"),
    make_option(c("-t", "--optimize_target"), type="character", default="ribo",
              help="optimize_target: ribo | te  [default= %default]", metavar="optimize_target"),
    make_option(c("-m", "--model_file"),type="character", default=NULL,
              help="path to the saved model file", metavar="model_file"),
    make_option(c("-j", "--n_cpus"), type="integer", default=16,
              help="number of cpus", metavar="n_cpus"),
    make_option(c("-o", "--out"), type="character", default="./output",
              help="output folder path [default= %default]", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

registerDoMC(cores=opt$n_cpus)


prefix=opt$prefix

annotFile=opt$annotation_file

model_fn=opt$model_file


feat.df=read.table(sprintf("%s.sparseFeature.txt.gz",prefix))
feat.mat=sparseMatrix(i=feat.df[,1],j=feat.df[,2],x=feat.df[,3],index1=F)
rownames(feat.mat)=read.table(sprintf("%s.sparseFeature.rowname",prefix),row.names=1)[,1]
colnames(feat.mat)=read.table(sprintf("%s.sparseFeature.colname",prefix),row.names=1)[1:ncol(feat.mat),1] ##because last few Feature(kmer) may not occurs in the data, so need to trim


TE.df=read.table(annotFile,row.names=1,header=T)


mRNA.RPKM.filter=TE.df[,'rpkm_rnaseq']>opt$min_rna_rpkm
ribo.RPKM.filter=TE.df[,'rpkm_riboseq']>opt$min_riboseq_rpkm
nomiss=complete.cases(TE.df)
TE.df=TE.df[ mRNA.RPKM.filter & nomiss,]



design_utr_len=100

print(load(model_fn))

ACGT=1:4
names(ACGT)=c("A","C","G","T")

###Kozak-EGFP GCCACC + 15bp cds
GFPcds="GCCACCATGGTGAGCAAGGGC"
flank1='TAAACTTAAGCTTGGTACCG'

featureExtraction<-function(seq){
cmd=sprintf("python  -W ignore FeatureExtraction_singleInput.py %s.sparseFeature.colname %s",prefix,paste(flank1,seq,GFPcds,sep=""))
data <- (read.table(pipe(cmd),sep=" ",header=F,comment.char=""))
return(unlist(data[1,]))
}

##should be produce logRibo level
if(opt$optimize_target=="ribo"){
predict_TE<-function(featVec){
predict(full.model.te,featVec)+predict(full.model.rna,featVec)
}
}else{
predict_TE<-function(featVec){
predict(full.model.te,featVec)
}

}



##ACGT => 1234
DNA2intVec<-function(seq){
unlist(lapply(unlist(strsplit(toupper(seq),"")),function(x) ACGT[x]))
}

##1234 => ACGT
intVec2DNA<-function(vec){
paste(unlist(lapply(vec,function(x) names(ACGT)[x])),collapse="")
}


######generate initial by high TE 5utr#####
fiveUTR_Population<-function(object){
popSize=object@popSize
cdslen=20
##use Claudia generated sequence as initial pool
raw.seqlist=read.table(file ="data/final_endogenous.txt")[,1]

len.list=unlist(lapply(raw.seqlist,function(x) nchar(x)-cdslen))
prob.list=len.list/sum(len.list)  ##give high chance to the longer UTR sequence
population.seq=sample(raw.seqlist, size=popSize, replace = T, prob = prob.list)
population=foreach(seq = population.seq,.combine=rbind)%do%{
seq=substr(seq,1,nchar(seq)-cdslen)
seq=gsub("N","",seq,ignore.case=T)
seq_design=substr(seq,21,nchar(seq))
DNA2intVec(seq_design)
}
 return(as.matrix(population))
}




Endogenous_maxFitness<-function(){
	raw.seqlist=read.table(file ="data/final_endogenous.txt")[,1]
	TE.list=foreach(seq=raw.seqlist)%dopar%{
	cmd=sprintf("python  -W ignore FeatureExtraction_singleInput.py %s.sparseFeature.colname %s",prefix,seq)
	data <- (read.table(pipe(cmd),sep=" ",header=F,comment.char=""))
	featVec=unlist(data[1,])
        predict_TE(featVec)
	}
	max(unlist(TE.list))
}

fitness<-function(intVec){
	utr.seq=intVec2DNA(intVec)
	featVec=featureExtraction(utr.seq)
	###avoid ATG in the UTR sequences
	if(grepl("ATG",utr.seq)||grepl("AT$", utr.seq)){  
		return(-10)
	}	
	predict_TE(featVec)
}


gaMutation_ACGT<-function(object, parent ){
    mutate <- parent <- as.vector(object@population[parent, ])
    n <- length(parent)
    j <- sample(1:n, size = 1)
    alphabet=1:4
    mutate[j] <- sample(alphabet[-mutate[j]],1)
    return(mutate)
}

randkey=sample(1:100000, 1)
logFile=sprintf("%s/%s.ga.%s.%d",opt$out,basename(model_fn),opt$optimize_target,randkey)
gaMonitor_writelog<-function(object, digits = getOption("digits")){

  fitness <- na.exclude(object@fitness)
    cat(paste("Iter =", object@iter, " | Mean =", format(mean(fitness),
        digits = digits), " | Best =", format(max(fitness), digits = digits),
        "\n"))
###write log file ###
##iter, seq, fitness
sink(logFile,append=T)
for(i in 1:nrow(object@population)){
	seq=intVec2DNA(object@population[i,])
	TE=object@fitness[i]
	cat(paste(object@iter,seq,TE,"\n"))
}
sink()
}

GA<-ga(type = "binary", fitness = fitness,nBits=design_utr_len ,lower = 1, upper = 4,maxiter=50, popSize=100,population=fiveUTR_Population,parallel=T,monitor=gaMonitor_writelog,mutation=gaMutation_ACGT,seed=randkey )



###Select diverge top Sequence from the evolution history ######


df=read.table(logFile)
colnames(df)<-c("iter",'seq','score')

###prune the similar sequence
nearestSeq<-function(seq1,topSeqlist){
	if(is.null(topSeqlist)){
		return(0)
	}
	scores=foreach(ts=topSeqlist)%do%{
		pairwiseAlignment(seq1, ts,scoreOnly=T)
	}
	which.max(unlist(scores))
}



prune2<-function(seq1,topSeqlist){
        if(is.null(topSeqlist)){
                return(100)
        }
        scores=foreach(ts=topSeqlist)%do%{
               # pairwiseAlignment(seq1, ts,scoreOnly=T,type="overlap")
        	unlist(adist(seq1,ts))
	}
        min(unlist(scores))
}



#####best N sequences must be better than the best endogenous sequence


similarThresh=160 #higher mean more similar
bestN=5
first3IterTop=1
first3IterBottom=1

selectBestNSequences<-function(){
topSeqlist=NULL
topSelID=NULL
bestInitScore=max(df[df[,'iter']==1,'score'])
df2=df[order(-df[,'score']),]
df2=df2[df2$score>bestInitScore,]
for(i in 1:nrow(df2)){
seq1=df2[i,'seq']

if(length(topSeqlist)>=bestN){
        break
}

if(prune2(seq1,topSeqlist)>5){
topSeqlist=c(topSeqlist,seq1)
topSelID=c(topSelID,i)
print(length(topSeqlist))
}

}

df2[unlist(topSelID),]
}

###bestN###
bestN.seqs=selectBestNSequences()

selectSimilarButBigChange<-function()
{
df1=df[df[,'iter']==1,] ##initial sequences
filter=df1[,'score']>-10
df1=df1[filter,]
df2=df[df[,'iter']==2|df[,'iter']==3,]
filter=df2[,'score']>-10
df2=df2[filter,]

scoreDlist=foreach(i=1:nrow(df2),.combine=rbind) %dopar%{
seq1=df2[i,'seq']
bestj=nearestSeq(seq1,df1[,'seq'])
scoreDiff=df2[i,'score']-df1[bestj,'score']
c(scoreDiff,bestj)
}
besti.top=which.max(scoreDlist[,1])
besti.bottom=which.min(scoreDlist[,1])
rbind(unlist(c(df2[besti.top,],sprintf("Up|%s|%f",df1[scoreDlist[besti.top,2],'seq'],df1[scoreDlist[besti.top,2],'score'] ))),
unlist(c(df2[besti.bottom,],sprintf("Down|%s|%f",df1[scoreDlist[besti.bottom,2],'seq'],df1[scoreDlist[besti.bottom,2],'score'] ))))
}

##similar but big change SBC
SBClist=selectSimilarButBigChange()
colnames(SBClist)[4]="description"

bestInitScore=max(df[df[,'iter']==1,'score'])

output.df=data.frame(bestN.seqs,description=rep(sprintf("bestTop|%f",bestInitScore),nrow(bestN.seqs)))
output.df=rbind(output.df,SBClist)
selSeqFile=sprintf("%s/%s.ga.%s.%d.sel",opt$out,basename(model_fn),opt$optimize_target,randkey)
write.table(output.df,file=selSeqFile,quote=F,col.names=F,sep="\t",row.names=F)
