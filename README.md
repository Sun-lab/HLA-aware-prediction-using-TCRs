# TCR-based Infection Prediction by Conditioning on HLA
This repository contains the source code associated with the paper, [Exploration of the Roles of HLAs When Predicting Infection Status by T Cell Receptors](https://www.biorxiv.org/content/10.1101/2024.11.18.624054v1).

A T cell can recognize some specific peptide-HLA complexes presented on cell surfaces by its T cell receptor (TCR). A blood sample can contain millions of unique TCRs, representing a snapshot of an individual’s TCR repertoire. TCR repertoire-wide association studies (TReWAS) aim to identify associations between individual TCRs and disease or exposure status. Prior research has shown that TCRs linked to viral infections can be discovered through TReWAS, and such TCRs can accurately predict current or past infections. Many of these TCRs are strongly associated with specific HLA alleles, suggesting that incorporating HLA information may enhance both TReWAS analyses and TCR-based predictions.

In our study, we evaluated TCR-based predictions while conditioning on individual HLA alleles or their k-nearest neighbors. We observed improved prediction accuracy for some HLA alleles. These HLA-specific predictions offer insights into the role of particular HLAs in infection or disease, with potential applications in personalized medicine.

# Terminology

A **HLA-agonistic Model** was trained using all individuals regardless of their HLA alleles. 

A **HLA-specific Model** for an HLA allele was trained using the individuals that carry a specific HLA allele. 

A **Combined Model** for an HLA allele was trained by combininig the TCRs identified by the HLA-agonistic model and the HLA-specific model.

A **KNN Model** for an HLA allele was trained by borrowing information from similar HLA alleles.

For details of how the models are trained and evaluated, please refer to the paper's method section.

# Requirements
There is no specific hardware requirements to excute all the scripts, Python scripts were developed under Python 3.11.13. For Python package versions please see environment.yml. 

# Directory Introduction

## Codes
The directory contains all the source codes to obtain the computational results, including codes to pre-process the data files (```codes/code_for_data```), codes for HLA specific model, HLA agnostic model, Combined model, and KNN model. The folder ```codes/HLA_CMV_association``` saves the source codes to compute the association result of each HLA with CMV. ```codes/intermediate_files``` saves the intermediate files needed during computations. ```codes/prepare_similarity_mat``` contains the files that prepares the similarity matrix which is needed for KNN method.

## Data
Required data files.

## Plots
The directory contains all the source code files to generate the plots in the paper. The sub-directory names correspond to each specific plots.

