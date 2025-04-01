# TCR-based Infection Prediction by Conditioning on HLA
This repository contains the source code associated with the paper, [Exploration of the Roles of HLAs When Predicting Infection Status by T Cell Receptors](https://www.biorxiv.org/content/10.1101/2024.11.18.624054v1).

A T cell can recognize some specific peptide-HLA complexes presented on cell surfaces by its T cell receptor (TCR). A blood sample can contain millions of unique TCRs, representing a snapshot of an individual’s TCR repertoire. TCR repertoire-wide association studies (TReWAS) aim to identify associations between individual TCRs and disease or exposure status. Prior research has shown that TCRs linked to viral infections can be discovered through TReWAS, and such TCRs can accurately predict current or past infections. Many of these TCRs are strongly associated with specific HLA alleles, suggesting that incorporating HLA information may enhance both TReWAS analyses and TCR-based predictions.

In our study, we evaluated TCR-based predictions while conditioning on individual HLA alleles or their k-nearest neighbors. We observed improved prediction accuracy for som HLA alleles. These HLA-specific predictions offer insights into the role of particular HLAs in infection or disease, with potential applications in personalized medicine.

# Terminology
A **HLA-agonistic Model** was trained using all individuals regardless of their HLA alleles. 

A **HLA-specific Model** for an HLA allele was trained using the individuals that carry a specific HLA allele. 

A **Combined Model** for an HLA allele was trained by combininig the TCRs identified by the HLA-agonistic model and the HLA-specific model.

For details of how the models are trained and evaluated, please refer to the paper's method section.

# Requirements
There is no specific hardware requirements to excute all the scripts, Python scripts were developed under Python 3.10.12. For Python package versions please see requirements.txt. 

# Directory Introduction

## Codes
The directory contains all the source codes to obtain the computational results, including codes to pre-process the data files (```code_for_data```), codes for ```HLA_specific_model```, ```HLA_agnostic_model```, and ```Combined_HLA_model```. The folder ```HLA_CMV_association``` saves the source dodes to compute the association result of each HLA with CMV. ```intermediate_files``` saves the intermediate files needed during computations. ```prepare_similarity_mat``` contains the files that prepares the similarity matrix which is needed for KNN method.

## Data
Required data files.

## Plots
The directory contains all the source code files to generate the plots in the paper. The sub-directory names correspond to each specific plots.

## Details for class II HLA allele similarity

The HLA class II allele similarity matrix is uploaded to github as the file:

[data/HLA_II_similarity_matrix_aa.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/HLA_II_similarity_matrix_aa.rds)

This file contains 250 HLA II alleles. What can be found in the HLA information for this group of 666 patients will be a subset of them. 
The order of corresponding HLAs in column is the same as that in the row names. 

Another detail that we mentioned in yesterday's meeting is about the haplotypes in the HLA information table:

[data/DeWitt_2018/HLA_v2_CMV_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/DeWitt_2018/HLA_v2_CMV_data.rds)

There are three types of more common HLA-II alleles, HLA-DRB1, HLA-DPAB, and HLA-DQAB. Each of these HLAs is decided by two chains, one alpha chain (A) and one beta chain (B), and that's how a name like ```HLA-DPAB*01:03_05:01``` comes from, which means that it is a pair ```(DPA*01:03, DPB*05:01)```. The names of HLA-DRB1 look different in format, for example,  ```HLA-DRB1*11:01```. This is due to that the alpha chain for HLA-DRA1 lacks diversity, and so only beta chain name is used to stand for the HLA allele. When doing analysis, each row in the HLA information file starting with "HLA-DRB1", "HLA-DPAB", and "HLA-DQAB" is treated as for one HLA allele, each decided by a pair of alpha chain and beta chain. The exception happens for 5 haplotypes, namely:

```
"HLA-DRDQ*03:01_05:01_02:01"
"HLA-DRDQ*09:01_03:02_03:03"
"HLA-DRDQ*10:01_01:05_05:01"
"HLA-DRDQ*13:01_01:03_06:03"
"HLA-DRDQ*15:01_01:02_06:02"
```

Each of these is a combination of one HLA-DRB1 and one HLA-DQAB allele. The authors of these data mentioned that this is due to, for example, in the first case, ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01``` always occur together among these individuals, and the authors combine them together to do the analysis. In our case, we can split each of them into two alleles, eg., split ```HLA-DRDQ*03:01_05:01_02:01``` into ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01```, and use the same corresponding occurrence information among the individuals for both of these two alleles. My understanding is that ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01``` have the exact the same occurrence among the individuals, but their pseudo amino acid sequences are different and thus will have different similarities with other HLA II alleles.  If we do not do clustering or k-nearest neighbor approach, the performances will be the same since they appear in the same group of individuals, but if we do clustering or k-nearest neighbor approach, the performances may be different since the pseudo amino acid sequence similarity is involved. 

For these haplotypes, our similarity matrix has data for each of the two components (eg., ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01``` ) but not for the haplotype format (eg., ```HLA-DRDQ*03:01_05:01_02:01```). 
