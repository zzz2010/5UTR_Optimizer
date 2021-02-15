#Author:Zhizhuo Zhang
#Date: 12/27/2015

.SECONDARY:
# output/gencode_v17_5utr_15bpcds.fa:extractSequence_5utr_cds.R
# 	Rscript $<

output/gencode_v17_5utr_15bpcds.fa.sparseFeature.txt.gz output/gencode_v17_5utr_15bpcds.fa.filter.fa:FeatureExtraction_final.py data/gencode_v17_5utr_15bpcds.fa
# 	use .viennarna-2.1.9
	python -W ignore $^


output/gencode_v17_5utr_15bpcds.fa.model:buildModel_final.R output/gencode_v17_5utr_15bpcds.fa.sparseFeature.txt.gz
	Rscript $< 

output/unique5UTR.predictedTE.txt:output/gencode_v17_5utr_15bpcds.fa.galog   output/gencode_v17_5utr_15bpcds.fa.claudia_seq.galog
	cat $^|tr ' '  '\t'|sort -k2|groupBy -g 2 -c 1,3 -o min,max|sort -k3gr > $@


output/unique5UTR.predictedTE%.txt:output/gencode_v17_5utr_15bpcds.fa%.galog   output/gencode_v17_5utr_15bpcds.fa%.claudia_seq.galog 
	cat $^|tr ' '  '\t'|sort -k2|groupBy -g 2 -c 1,3 -o min,max|sort -k3gr > $@



all_evojob.TE:
	seq 1 100|Rscript evolutionDesign_TE.R muscle
	seq 1 100|Rscript evolutionDesign_TE.R pc3 
	seq 1 100|Rscript evolutionDesign_TE.R HEK_Andrev2015 
	
	
all_evojob.Ribo:	
	seq 1 100|Rscript evolutionDesign_Ribo.R muscle 
	seq 1 100|Rscript evolutionDesign_Ribo.R pc3 
	seq 1 100|Rscript evolutionDesign_Ribo.R HEK_Andrev2015 
 
	
all_seljob:
	find output/final/raw/|grep fa|xargs -n 1 -I {}  Rscript select_diverseDesginSequence.R {} 
 

output/final/synthetic3K.txt:finalFormat_3k_synthetic_seqs.py
	python $< > $@
