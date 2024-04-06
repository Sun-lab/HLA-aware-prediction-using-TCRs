# Conditional TCR Prediction

Use TCR to predict infection history, while conditioning on HLA.

## To-dos

1. Illustration that borrow information could help. We can just use the aa distance for now. 

    - Conduct herarchical clustering of HLAs. For each HLA, we have three sets of results. (1) the results of HLA-specific prediction. (2) The results for this HLA using the model trained for the HLA cluster to which this HLA belongs. (3) The results for this HLA using the model trained by all HLAs. For (2), we just evaluate the model using part of testing data: the individuals who have the HLA of interest. For (3), we can train 20 models for all HLAs with 20 random splits of training/testing data. Then when we want to evaluate one HLA, we can use the individuals in the testing data who have the HLA of interest. Therefore, we do not need to retrain the model of all HLAs for each HLA cluster. 

2. Instead of clustering, we can find the neighbors for each HLA, using a distance cutoff. Then we try to make prediction by borrowing information across HLAs. We may give different weights for different HLAs. For example, if we are intersted in HLA_1, and we got two neighbors: HLA_2 and HLA_3, with similarity 0.8 and 0.6. Then we may take the sampels with any of the three HLAs, but give them diffrent weights. For example, weight 1 for HLA_1, 0.8^2 and 0.6^2 for the two other HLA. We can think about what is the optimal weights here, or what is the optimal transformation weight = f(similarity) where f should be a monontone function.

   - We can use weights to find TCR-CMV associations by a logistic regression, where can supply the weight for each observation. When using logistic regression, we need to use LRT to obtain p-value because the data are highly sparse and the default Wald test has inflated type I error.
   - If we define y = 1 or 0 for the appearance of a TCR, and x = 1 or 0 for the case/control status (CMV or not), rd is the total number of TCRs per sample. 
   - model 0: y ~ log(rd)
   - model 1: y ~ log(rd) + x
   - anova(m0, m1, test = "LRT").
  
3. Examine some examples where clustering has much better results than weighted or the other way around.
   - Check whether the results become consistent when we use the same train/test split, 500 training and 164 for testing.
   - Compare resutls for each HLA for 4 methods. model trained using all individuals, only individual with this HLA, clustering and k-nearest neighbors. 

   
## Finished

1. Collect TCR data, HLA data, and CMV infection status.

   TCR data: [data/Emerson_2017/TCR_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/Emerson_2017/TCR_data.rds)

   HLA data and CMV data: (the row at the bottom is CMV. NA is missing value.)[data/DeWitt_2018/HLA_v2_CMV_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/DeWitt_2018/HLA_v2_CMV_data.rds)

   In both files, each column is one of the 666 individual, and the order of the columns is consistent cross the two files.

2. We want to use the total numbers of CMV-associted TCRs to predict CMV status. Evaluate the prediction using AUC in testing data. We have the data of y_i = CMV status of individual i, x_i = # of TCRs of individual i that are associated CMV at a particular p-value cutoff, x_i/d_i, where d_i is the total number of TCRs of individual i.

3. For conditional analysis. Just use those individuals with certain HLA to conduct the above analysis.


    Data file of HLA similarity matrices:

    [data/correlation_matrices.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/correlation_matrices.rds)

    This file is a list of four matrices, and we may care about three of them for now. Their names have correspondence as explained below:

    "aa" -> "AA" in Fi. 6(B) in the grant -> similarity computed purely based on amino acids in HLA pseudo sequences

    "cmv_10000" -> "CMV" in Fig. 6(A) in the grant -> similarity computed as the correlation of prediction cores given by DePTH model (a model predicting the association of TCR and HLA given their amino acid sequences) between HLA and top CMV associated TCRs

    "covid_29593" -> SARS-COV-2 in Fig. 6(A, B) in the grant -> similarity computed as the correlation of prediction cores given by DePTH model (a model predicting the association of TCR and HLA given their amino acid sequences) between HLA and covid associated TCRs

    For each matrix, the row names give the corresponding HLA. The column names are not specified, but they follow the same order as the row names.

    One way of getting a distance matrix is to do 1-similarity matrix:
   ```dist1 = as.dist(1 - similarity_matrix); h1 = hclust(d=dist1)```



Try to reproduce the Figure 7 in the grant.

## Details for class II HLA allele similarity

The HLA class II allele similarity matrix is uploaded to github as the file:

[data/HLA_II_similarity_matrix_aa.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/HLA_II_similarity_matrix_aa.rds)

This file contains 250 HLA II alleles. What can be found in the HLA information for this group of 666 patients will be a subset of them. 
The order of corresponding HLAs in column is the same as that in the row names. 

Another detail that we mentioned in yesterday's meeting is about the haplotypes in the HLA information table:

[data/DeWitt_2018/HLA_v2_CMV_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/DeWitt_2018/HLA_v2_CMV_data.rds)

There are three types of more common HLA-II alleles, HLA-DRB1, HLA-DPAB, and HLA-DQAB. Each of these HLAs is decided by two chains, one alpha chain (A) and one beta (chain), and that's how a name like ```HLA-DPAB*01:03_05:01``` comes from, which means that it is a pair ```(DPA*01:03, DPB*05:01)```. The names of HLA-DRB1 look different in format, for example,  ```HLA-DRB1*11:01```. This is due to that the alpha chain for HLA-DRA1 lacks diversity, and so only beta chain name is used to stand for the HLA allele. When doing analysis, each row in the HLA information file starting with "HLA-DRB1", "HLA-DPAB", and "HLA-DQAB" is treated as for one HLA allele, each decided by a pair of alpha chain and beta chain. The exception happens for 5 haplotypes, namely:

```
"HLA-DRDQ*03:01_05:01_02:01"
"HLA-DRDQ*09:01_03:02_03:03"
"HLA-DRDQ*10:01_01:05_05:01"
"HLA-DRDQ*13:01_01:03_06:03"
"HLA-DRDQ*15:01_01:02_06:02"
```

Each of these is a combination of one HLA-DRB1 and one HLA-DQAB allele. The authors of these data mentioned that this is due to, for example, in the first case, ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01``` always occur together among these individuals, and the authors combine them together to do the analysis. In our case, we can split each of them into two alleles, eg., split ```HLA-DRDQ*03:01_05:01_02:01``` into ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01```, and use the occurrence information among the individuals for both of these two alleles. My understanding is that ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01``` have the exact the same occurrence among the individuals, but their amino acid sequences are different and thus will have different similarities with other HLA II alleles.  If we do not do clustering or k-nearest neighbor approach, the performances will be the same since they appear in the same group of individuals, but if we do clustering or k-nearest neighbor approach, the performances may be different since the amino acid sequence similarity is involved. 

For these haplotypes, our similarity matrix has data for each of the two components (eg., ```HLA-DRB1*03:01``` and ```HLA-DQAB*05:01_02:01``` ) but not for the haplotype format (eg., ```HLA-DRDQ*03:01_05:01_02:01```). 
