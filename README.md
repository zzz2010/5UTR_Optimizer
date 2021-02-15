# Design of evolved synthetic 5'UTR sequences for enhanced protein production

The goal of this workflow is to generate novel 5'UTR sequences with optimized features to maximize the translation efficiency of the constructs. 
This is achieved in 3 steps:

- Step 1. **Feature extraction** from endogenous sequences 
- Step 2. **Random Forest model generation** and evaluation, using extracted features to predict Translation Efficiency (TE)
- Step 3. Generation of novel **'evolved' sequences using a Genetic Algorithm** 

This repository is the code accompanying the paper:High-Throughput 5’ UTR Engineering for Enhanced Protein Production in Non-Viral Gene Therapies. Jicong Cao*, Eva Maria Novoa*, Zhizhuo Zhang*, William C.W. Chen, Dianbo Liu, Gigi C G Choi, Alan S L Wong, Claudia Wehrspaun, Manolis Kellis, Timothy K Lu. bioRxiv 2020. doi: https://doi.org/10.1101/2020.03.24.006486


![alt text](./img/init_fig.png "init_fig")


## Table of Contents
- [Step 1. Design of 5'UTR sequences](#Step-1.-Design-of-5'UTR-sequences)
- [Step 2. 5'UTR model generation and evaluation](#Step-2.-5'UTR-model-generation-and-evaluation)
- [Step 3. Generation of novel evolved 5'UTR sequences](#Step-3.-Generation-of-novel-evolved-5'UTR-sequences)
- [Running individual scripts for each of the steps performed](#Running-individual-scripts-for-each-of-the-steps-performed)
- [Additional data sources](#Additional-data-sources)
- [Dependencies and versions](#Dependencies-and-versions)
- [Citation](#Citation) 
- [Contact](#Contact) 
 

## Step 1. Feature extraction of 5'UTR sequences 

**Goal**: Extract features that correspond to optimal 5'UTRs to maximize the translation efficiency (TE).

Features used in Random Forest model training:
- Kmer-frequency: k=1-6
- RNAfolding energy for first 100bp, last 30bp(15bp UTR+ 15bp CDS), whole (5UTR+15bp CDS), 5UTR only, with/without consdiering G-quadruplex 
- codon usage
- number of start and stop codon in UTR
- 5'UTR length

### Running the code: 
``` 
python run_pipeline.py feature_extract --input_fasta data/gencode_v17_5utr_15bpcds.fa --output_dir output/
```


## Step 2. 5'UTR model generation and evaluation

**Goal**: Generate and evaluate model  features that correspond to optimal 5'UTRs to maximize the translation efficiency (TE).

**Algorithm used**: **Random Forest** to build a prediction model 

The model trained for human sequences is available as part of this repository. 
Results: Using 10-fold cross-validation, we obtained  0.71 pearson correlation in TE prediction, and 0.74 in RNA expression prediction.

### Running the code: 
Build full model
``` 
python run_pipeline.py model_build --prefix output/input.fa --annotation_file data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt --min_rna_rpkm 5 --min_riboseq_rpkm 0.1 --model 1 --out output/muscle_randomforest.model
```
Evaluate different choices of ML model
``` 
python run_pipeline.py model_eval --prefix output/input.fa --annotation_file data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt --min_rna_rpkm 5 --min_riboseq_rpkm 0.1 --model 1 --out output/muscle_randomforest.model
```

## Step 3. Generation of novel evolved 5'UTR sequences

**Goal**: Evolve endogenous 5'UTR sequences to obtain 5'UTRs with increased translation efficiency

**Algorithm used**: **Genetic Algorithm (GA)**

Additional details of individual steps are described in `makefile` file)

#### Running the code:
```
python run_pipeline.py  sequence_generate --n_total 3585 --prefix output/input.fa --annotation_file data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt --min_rna_rpkm 5 --min_riboseq_rpkm 0.1 -t ribo  -m output/muscle_randomforest.model  -o output/
```
 


## Additional data sources 

- df_counts_and_len.TE_sorted.with_annot.txt: Muscle RNA-seq and Ribo-seq data from publicly available datasets (PC3, HEK and muscle)
- entiredata_2015_12_25_3-25_pm: RNA binding motif, downloaded from http://cisbp-rna.ccbr.utoronto.ca/
 
## Dependencies and versions

#### Create a Conda Environment with all the dependencies:
```
conda env create -f environment.yml
```

Here we list the specific versions of individual dependencies needed:

| Software | Version |
| ------------- | ------------- |
| python  | 3.5  |
| viennarna  | 2.1.9  |
| R  | 3.5.1  |
| biopython  | 1.72 |
| r-randomforest  | 4.6_14 |
| r-ga  | 3.2  |
| r-seqinr  | 3.6_1  |
| r-glmnet  | 2.0_16  |


## Citation
If you find this work useful, please cite: 

High-Throughput 5’ UTR Engineering for Enhanced Protein Production in Non-Viral Gene Therapies. Jicong Cao*, Eva Maria Novoa*, Zhizhuo Zhang*, William C.W. Chen, Dianbo Liu, Gigi C G Choi, Alan S L Wong, Claudia Wehrspaun, Manolis Kellis, Timothy K Lu. bioRxiv 2020. doi: https://doi.org/10.1101/2020.03.24.006486


## Contact
If you have any issues using this code, please open an Issue in the GitHub repository. Thanks!


#### Last Updated: February 2021

