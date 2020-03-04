###Author: Zhizhuo Zhang
###Date: 12/27/2015
###use bigger training data 
###smaller pop size, and fewer generations
###run more rounds


library(GA)
library(randomForest)
library(methods)
library(Matrix)
library(foreach)
library(doMC)
library(Metrics)
library(seqinr)
options(stringsAsFactors = FALSE)

registerDoMC(cores=30)

prefix="output/gencode_v17_5utr_15bpcds.fa"

cell=".HEK_Andrev2015" ##HEK cell
cell=".pc3"
#cell=""  ## muscle 

Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]
cell2=Args[2]
if(cell2=="muscle"){
cell=""
}else{
cell=paste(".",cell2,sep="")
}

TE.df=read.table(sprintf("data/df_counts_and_len.TE_sorted%s.with_annot.txt",cell),row.names=1,header=T)
if(cell==""){
mRNA.RPKM.filter=TE.df[,'rpkm_rnaseq']>5
ribo.RPKM.filter=TE.df[,'rpkm_riboseq']>0.1
}else{
mRNA.RPKM.filter=TE.df[,'rpkm_rnaseq']>50
ribo.RPKM.filter=TE.df[,'rpkm_riboseq']>5
}


nomiss=complete.cases(TE.df)
TE.df=TE.df[ mRNA.RPKM.filter & nomiss,]



design_utr_len=100

if(cell!=""){
cell=paste(cell,".big",sep="")
}
print(cell)
print(load(sprintf("%s%s.model",prefix,cell)))

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
predict_TE<-function(featVec){
predict(full.model.te,featVec)+predict(full.model.rna,featVec)
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
raw.seqlist=read.table(file ="output/final_endogenous.txt")[,1]

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
	raw.seqlist=read.table(file ="output/final_endogenous.txt")[,1]
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

gaMonitor_writelog<-function(object, digits = getOption("digits")){
logFile=sprintf("output/final/raw/%s.%s.final.gaRibo.%d",basename(prefix),cell2,randkey)
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

GA<-ga(type = "binary", fitness = fitness,nBits=design_utr_len ,min = 1, max = 4,maxiter=50, popSize=100,population=fiveUTR_Population,parallel=T,monitor=gaMonitor_writelog,mutation=gaMutation_ACGT,seed=randkey )

