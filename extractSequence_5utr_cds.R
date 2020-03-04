###Author: Zhizhuo Zhang
###Date: 12/27/2015
library(parallel)
library(doMC)
registerDoMC(cores=30)
library("BSgenome.Hsapiens.UCSC.hg19") # here use whatever annotation of version of genome instead - eg. claudia is using hg38, but Eva using hg19 since the reads were mapped in hg19
genome<-BSgenome.Hsapiens.UCSC.hg19 
library(GenomicFeatures)
geneCacheFn="~/compbio/share_data/genomes/hg19db_GencodeV17.sqlite"

txdb <- loadDb(geneCacheFn)
cds<-unlist(cdsBy(txdb, by="tx", use.names=TRUE))
fiveUTRs<-unlist(fiveUTRsByTranscript(txdb,use.names=TRUE))				

length(fiveUTRs)
###filter very short 5UTR ###
#filter=unlist(mclapply(fiveUTRs,function(x) width(ranges(x))>15))
#fiveUTRs=fiveUTRs[filter]
#length(fiveUTRs)

#get the first cds
match.cds=nearest(fiveUTRs,cds,select="all")
library(foreach)
topN=length(match.cds)
bed12format=foreach(k=1:topN,.combine=rbind)%dopar%{
i=match.cds[k]@queryHits
j=match.cds[k]@subjectHits
tx.id1=names(fiveUTRs[i])
tx.id2=names(cds[j])

#make sure they are the same transcript
if(tx.id1!=tx.id2){
return(NULL)
}

strand=as.character(strand(fiveUTRs[i]))

utr.r=ranges(fiveUTRs[i])
cds.r=ranges(cds[j])
##take first 15bp of CDS based on strand##
##order bed12 exon list based on strand ##
if(width(utr.r)+width(cds.r)<30){
	return(NULL)
}

chr=as.character(seqnames(fiveUTRs[i]))
start=min(start(utr.r), start(cds.r))
end=max(end(utr.r),end(cds.r))

firstNbase=15
if(strand=="+"){
	width(cds.r)=firstNbase
	blocksize.str=sprintf("%d,%d",width(utr.r) ,width(cds.r)  )
	blockstart.str=sprintf("%d,%d",utr.r@start-start, cds.r@start-start )
}else{
	start(cds.r)=start(cds.r)+width(cds.r)-firstNbase
	blocksize.str=sprintf("%d,%d",width(cds.r),width(utr.r))
	blockstart.str=sprintf("%d,%d",cds.r@start-start , utr.r@start-start)
}


#name=sprintf("%s__%d_%d",tx.id1,start,end)
c(chr,start,end,tx.id1,1,strand,start,end,"0,0,0",2,blocksize.str, blockstart.str  )

}


write.table(bed12format,"output/5utr_1stcds.bed12",sep="\t",col.names=F,quote=F,row.names=F)

cmd="fastaFromBed -fi ~/compbio/share_data/genomes/hg19.fa -bed output/5utr_1stcds.bed12  -name -s -split -fo output/gencode_v17_5utr_15bpcds.fa"
system(cmd)
