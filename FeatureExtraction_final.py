####extract feature for latter prediction use
###RNA folding
###Condon
### k-mer
### motif
import sys,os
from Bio import SeqIO
import Bio.SeqUtils.CodonUsage
import subprocess
from multiprocessing import Pool,cpu_count
import gzip
from FeatureCommons import *

cds_length=15 ##assumption last portion of sequence is cds

#inputFasta="output/test.fa"
inputFasta=sys.argv[1]
out_dir=sys.argv[2] ##output directory


###take the longest length for the same transcript id
tx_seq=dict()
for seq_record in SeqIO.parse(inputFasta, "fasta"):
    tx=seq_record.id
    seq=seq_record.seq
    if len(seq)<30: #skip when it too short
        continue
    if "ATG" not in seq:
        continue
    if tx not in tx_seq or len(seq)>len(tx_seq[tx]):
        tx_seq[tx]=seq
    

###output the non-redundancy fasta
outputFasta=out_dir+"/input.filter.fa"
outf=open(outputFasta,"w")
txIDlist=list()
seqList=list()
for tx in tx_seq:
    outf.write(">"+tx+"\n")
    outf.write(tx_seq[tx].tostring()+"\n")
    txIDlist.append(tx)
    seqList.append(tx_seq[tx])
outf.close()




pool = Pool(cpu_count())
featList=pool.map(Seq2Feature,seqList)

##output feature matrix, tx.id, feature.id, feature.value, id is 0-based
outf2=gzip.open(out_dir+"/input.fa.sparseFeature.txt.gz",'wt')
feat2ID=dict()
featid=-1
for i in range(len(txIDlist)):
    txid=i
    for featItem in featList[i]:
        featname=featItem[0]
        featVal=featItem[1]        
        if featname not in feat2ID:
            featid+=1
            feat2ID[featname]=featid
            fid=featid
        else:
            fid=feat2ID[featname]
        outstr=str(i)+"\t"+str(fid)+"\t"+str(featVal)
        outf2.write(outstr+"\n")

outf2.close()

##mapping id to human understandable name
outf3=open(out_dir+"/input.fa.sparseFeature.rowname",'w')
for i in range(len(txIDlist)):
    outf3.write(str(i)+"\t"+txIDlist[i]+"\n")

outf3.close()

outf4=open(out_dir+"/input.fa.sparseFeature.colname",'w')
sorted_items=sorted(feat2ID.items(), key=lambda x: x[1])
for a in sorted_items:
    print(a)
    featname=a[0]
    fid=a[1]
    outf4.write(str(fid)+"\t"+featname+"\n")
outf4.close()



