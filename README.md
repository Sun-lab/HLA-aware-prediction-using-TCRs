# Conditional TCR Prediction

Use TCR to predict infection history, while conditioning on HLA.

## To-dos

1. Collect TCR data, HLA data, and CMV infection status.

   TCR data (TCRs are already filtered. The kept TCRs should all each appears in at least 7 among the all 666 individuals. There should not be any missing value. ):
   HLA data and CMV data (the rows on the top corresponds to HLAs and the row at the end is CMV. NA is missing value.):

2. Split the data to training and testing, maybe 4 folds for training and 1 fold for testing.
   2.1 Just sample 80% of individuals as training samples. sample(666, round(0.8*666)).

3. In the training data
   3.1 Assess the association between TCRs and CMV status, using Fisher's exact test. We may only use the TCRs appear in at least 5 or 7 individuals.

4. Select the TCRs with p-values (from training data) smaller than a threshold, and use the total numbers of such TCRs to predict CMV status in both training and test data. Evaluate the prediction using AUC. We have the data of y_i = CMV status of individual i, x_i = # of TCRs of individual i that are associated CMV at a particular p-value cutoff, x_i/d_i, where d_i is the total number of TCRs of individual i.

5. For conditional analysis. Just use those individuals with certain HLA to conduct the above analysis.

Try to reproduce the Figure 7 in the grant.
