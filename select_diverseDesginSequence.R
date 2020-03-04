
similarThresh=160 #higher mean more similar
bestN=5
first3IterTop=1
first3IterBottom=1
library("seqinr")
library(Biostrings)
library(foreach)
library(doMC)
registerDoMC(cores=30)

Args<-commandArgs()[grep("^--",commandArgs(),invert=T)]
inputFn="output/final/raw/gencode_v17_5utr_15bpcds.fa.muscle.claudia_seq.galog.18870"
inputFn=Args[2]

options(stringsAsFactors = FALSE)

df=read.table(inputFn)
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

write.table(output.df,file=sprintf("output/final/sel/%s",basename(inputFn)),quote=F,col.names=F,sep="\t",row.names=F)
