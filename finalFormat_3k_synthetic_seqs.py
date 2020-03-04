import os,sys
import glob
import math


allFiles=glob.glob("output/final/sel/*")


model_ostr_score=dict()
visited=set()
#gencode_v17_5utr_15bpcds.fa.pc3.galog.62355
#gencode_v17_5utr_15bpcds.fa.pc3.claudia_seq.gaRibo.13167#
Nprinted=0
print "seq	score	model	generation	info"
for fn in allFiles:
	label=os.path.basename(fn).replace("gencode_v17_5utr_15bpcds.fa.","").replace(".claudia_seq","").replace(".gaRibo","_Ribo").replace(".galog","_TE").replace(".final","").split(".")[0]
	for line in open(fn):
		comps=line.strip().split()
		seq="TAAACTTAAGCTTGGTACCG"+comps[1]+"GCCACCATGGTGAGCAAGGG"
		if seq in visited:
			continue
		visited.add(seq)
		score=comps[2]
		if score=="NA":
			continue
		itera=comps[0]
		info=comps[3]
		outstr=seq+"\t"+score+"\t"+label+"\t"+itera+"\t"+info
		if info.startswith("best"):
			initBestScore=float(info.split("|")[1])
			if float(score)<initBestScore+0.05:
				continue
			if label not in model_ostr_score:
				model_ostr_score[label]=dict()
			model_ostr_score[label][outstr]=float(score)
		else:
			print outstr
			Nprinted+=1
Ntotal=3585
import operator
perModelNum=int(math.ceil(float(Ntotal-Nprinted)/len(model_ostr_score)))

for model in model_ostr_score:
	x=model_ostr_score[model]
	sorted_x = sorted(x.items(), key=operator.itemgetter(1),reverse=True)
	for i in range(perModelNum):
		print sorted_x[i][0]
