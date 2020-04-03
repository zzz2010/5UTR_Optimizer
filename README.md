# Code accompanying paper:

Jicong Cao, Eva Maria Novoa, Zhizhuo Zhang, William C.W. Chen, Dianbo Liu, Gigi C G Choi, Alan S L Wong, Claudia Wehrspaun, Manolis Kellis, Timothy K Lu. High-Throughput 5’ UTR Engineering for Enhanced Protein Production in Non-Viral Gene Therapies. bioRxiv 2020. doi: https://doi.org/10.1101/2020.03.24.006486


1.Design of 5'UTR sequences
==========

Goal: Design optimal 5'UTRs to maximize the translation efficiency (TE).
Algorithm used: **Random Forest** to build a prediction model 

Features used in Random Forest model training:
- Kmer-frequency: k=1-6
- RNAfolding energy for first 100bp, last 30bp(15bp UTR+ 15bp CDS), whole (5UTR+15bp CDS), 5UTR only, with/without consdiering G-quadruplex 
- codon usage
- number of start and stop condon in UTR
- 5UTR length

The model seems to work well. Using 10fold cross-valiation, we obtained  0.71 pearson correlation in TE prediction, and 0.74 in RNA expression prediction

2.Workflow (detail is described in `makefile` file)
=====================

Goal: evolve endogenous 5'UTR sequences to obtain 5'UTRs with increased translation efficiency
Algorithm used: Genetic Algorithm (GA)

Steps performed: 
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


3.Data sources
===========
- df_counts_and_len.TE_sorted.with_annot.txt: Muscle RNA-seq and Ribo-seq data from publicly available datasets (PC3, HEK and muscle)
- entiredata_2015_12_25_3-25_pm: RNA binding motif, downloaded from http://cisbp-rna.ccbr.utoronto.ca/
 

### Last Updated: April 2020
