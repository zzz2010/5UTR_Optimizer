####extract feature for latter prediction use
###RNA folding
###Condon
### k-mer
### motif
import sys,os
from Bio import SeqIO
import Bio.SeqUtils.CodonUsage
import subprocess
from multiprocessing import Pool
import gzip
from Bio.Seq import Seq
from FeatureCommons import *
cds_length=15 ##assumption last portion of sequence is cds

#featureDefineFn="output/gencode_v17_5utr_15bpcds.fa.sparseFeature.colname" #sys.argv[1]
#seq="TAGCCCCAAATCCTTGGGCGGCTGGGAGGTGGTCCCTGTCATGCCCACGACGGGTGGGCGCGTCTTTACCGACAGCCTGCATGGATCCGCAGATGTGCTGGGCGCTTAGATGCCTGCGGCTATGCAGAGACGCAGTCCAGAGAGGTGGGTAAGATCAGGGGA" #sys.argv[2]

featureDefineFn=sys.argv[1]
seq=sys.argv[2]
seq=Seq(seq)

feature_id=dict()
for line in open(featureDefineFn):
	comps=line.strip().split()
	feature_id[comps[1]]=int(comps[0])


featList=Seq2Feature(seq)

FeatureVec=[0]*len(feature_id)
for item in featList:
	i=feature_id[item[0]]
	val=item[1]
	FeatureVec[i]=val

print(" ".join(map(str,FeatureVec)))


