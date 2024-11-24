# Conditional TCR Prediction

Use TCR to predict infection history, while conditioning on HLA.

## Codes
The directory contains all the source code files to obatin the computational results.

### Codes/Combined_HLA_model
This directoy contains the source code file to train and make predictions for the Combined Model in the paper.

### Codes/HLA_CMV_association
This directoy contains the source code file to compute the association result of each HLA with CMV. 

### Codes/HLA_agnostic_model



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
