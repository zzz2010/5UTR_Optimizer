### Author:Zhizhuo Zhang
### Date: 12/27/2015


5utrdesign
==========

Design optimal 5utr to maximize the translation efficiency.
I use Randomforest to build a prediction model with the following features:
- Kmer-frequency: k=1-6
- RNAfolding energy for first 100bp, last 30bp(15bp UTR+ 15bp CDS), whole (5UTR+15bp CDS), 5UTR only, with/without consdiering G-quadruplex 
- codon usage
- number of start and stop condon in UTR
- 5UTR length

The model seems work well. Using 10fold cross-valiation, I got 0.71 pearson correlation in TE prediction, and 0.74 in RNA expression prediction

Workflow (detail is described in `makefile` file)
=====================
- extract DNA sequence 5'UTR+first CDS
> make output/gencode_v17_5utr_15bpcds.fa

- compute sequence feature
> use .viennarna-2.1.9
> use .biopython-1.64-python-2.7.1-sqlite3-rtrees
> make output/gencode_v17_5utr_15bpcds.fa.sparseFeature.txt.gz

- build prediction model for TE and Ribo-seq expression
> make output/gencode_v17_5utr_15bpcds.fa.model

- design optimal 100bp 5UTR sequence for maximizing TE
> make all_evojob.TE

- design optimal 100bp 5UTR sequence for maximizing Ribo-seq expression
> make all_evojob.Ribo

- select diverse optimized sequences for UTR synthesis 
> make all_seljob

> make output/final/synthetic3K.txt


Data
===========
- df_counts_and_len.TE_sorted.with_annot.txt: Muscle RNA-seq and Ribo-seq data from Eva Maria Novoa
- entiredata_2015_12_25_3-25_pm: RNA binding motif, downloaded from http://cisbp-rna.ccbr.utoronto.ca/
 
