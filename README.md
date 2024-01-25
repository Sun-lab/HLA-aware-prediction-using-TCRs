# Conditional TCR Prediction

Use TCR to predict infection history, while conditioning on HLA.

## To-dos

1. Collect TCR data, HLA data, and CMV infection status.

   TCR data (TCRs are already filtered. The kept TCRs should all each appears in at least 7 among the all 666 individuals. There should not be any missing value. ):

   [data/Emerson_2017/TCR_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/Emerson_2017/TCR_data.rds)

   HLA data and CMV data (the rows on the top corresponds to HLAs and the row at the bottom is CMV. NA is missing value.):

   [data/DeWitt_2018/HLA_v2_CMV_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/DeWitt_2018/HLA_v2_CMV_data.rds)

   In both files, each column is one of the 666 individual, and the order of the columns is consistent cross the two files.

2. Split the data to training and testing, maybe 4 folds for training and 1 fold for testing.
   2.1 Just sample 80% of individuals as training samples. sample(666, round(0.8*666)).

3. In the training data
   3.1 Assess the association between TCRs and CMV status, using Fisher's exact test. We may only use the TCRs appear in at least 5 or 7 individuals.

4. Select the TCRs with p-values (from training data) smaller than a threshold, and use the total numbers of such TCRs to predict CMV status in both training and test data. Evaluate the prediction using AUC. We have the data of y_i = CMV status of individual i, x_i = # of TCRs of individual i that are associated CMV at a particular p-value cutoff, x_i/d_i, where d_i is the total number of TCRs of individual i.

5. For conditional analysis. Just use those individuals with certain HLA to conduct the above analysis.
   - 5.1 We need to repeat the split of training/testing for each HLA multiple times, to record the mean and sd of AUC. 

6. Try to borrow information across HLAs to improve the prediction accuracy for a specific HLA.

    - 6.1 Conduct herarchical clustering of HLAs using aa distance of cmv_10000 distance. Cut the tree at certain distance to select small clusters of HLAs.
    - 6.2 For each cluster, try to take the union of the samples who have any of those HLA alleles in the clsuter, and make prediction for them.
    - 6.3 To compare the results of HLA clusters, we can first compare with the results of individual HLA (using the results of 5.1). We should also compare with the results using all HLAs. In this case, we can just train 20 models for all HLAs. Then when we want to evaluate one cluster, we can use part of the testing data that are individuals who have any HLA in the cluster. Therefore, we do not need to retrain the model of all HLAs for each HLA cluster. 
    - 6.4 Repeat the analysis for aa distance. 

7. When we want to improve the prediction for one HLA, we may find a few neighboring HLA to help. When we build the model, we may give different weights for different HLAs. For example, if we are intersted in HLA_1, and we got two neighbors: HLA_2 and HLA_3, with similarity 0.8 and 0.6. Then we may take the sampels with any of the three HLAs, but give them diffrent weights. For example, weight 1 for HLA_1, 0.8^2 and 0.6^2 for the two other HLA. We can think about what is the optimal weights here, or what is the optimal transformation weight = f(similarity) where f should be a monontone function. 

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
