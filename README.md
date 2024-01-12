# Conditional TCR Prediction

Use TCR to predict infection history, while conditioning on HLA.

## To-dos

1. Collect TCR data, HLA data, and CMV infection status.

   TCR data (TCRs are already filtered. The kept TCRs should all each appears in at least 7 among the all 666 individuals. There should not be any missing value. ):
   
   [data/Emerson_2017/TCR_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/Emerson_2017/TCR_data.rds)
   
   HLA data and CMV data (the rows on the top corresponds to HLAs and the row at the bottom is CMV. NA is missing value.):
   
   [data/DeWitt_2018/HLA_v2_CMV_data.rds](https://github.com/Sun-lab/conditional_TCR_prediction/blob/main/data/DeWitt_2018/HLA_v2_CMV_data.rds)

   In both files, each column is one of the 666 individual, and the order of the columns is consistent cross the two files.
   
3. Split the data to training and testing, maybe 4 folds for training and 1 fold for testing.
   2.1 Just sample 80% of individuals as training samples. sample(666, round(0.8*666)).

4. In the training data
   3.1 Assess the association between TCRs and CMV status, using Fisher's exact test. We may only use the TCRs appear in at least 5 or 7 individuals.

5. Select the TCRs with p-values (from training data) smaller than a threshold, and use the total numbers of such TCRs to predict CMV status in both training and test data. Evaluate the prediction using AUC. We have the data of y_i = CMV status of individual i, x_i = # of TCRs of individual i that are associated CMV at a particular p-value cutoff, x_i/d_i, where d_i is the total number of TCRs of individual i.

6. For conditional analysis. Just use those individuals with certain HLA to conduct the above analysis.

7. Try to borrow information across HLAs to improve the prediction accuracy for a specific HLA.
  7.1 Calculate HLA similarity using their pseudo sequence using Blossum62. (Si to share the file if this has been calculated).
  7.2 Calculate HlA similarity using their associations with TCRs. Need to first predict HLA-TCR associations. Then the data for one HLA is its assocaition with, say 10,000 TCRs. 

Try to reproduce the Figure 7 in the grant.
