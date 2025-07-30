specific_HLA_pred, specific_HLA_II_pred

	split the 85 HLAs to 9 chunks, run them by parallele jobs, sample_ids.txt

	specific_HLA_pred.py is the main code. 


all_sample_pred

	HLA agnostic model

	all_sample_pred.py: the main code to run the prediction using all samples. 

	diff_pvals.py: use different p-value cutoffs.
		constrained_pred make prediction only for the individuals with a specific HLA. 


all_model_split
	
	For HLA with sample size between 40 and 70, we spllit the data multiple times and try to average or aggregate the results

weighted_pred, weighted_pred_HLAII
	
	KNN prediction when we weight the HLAs based on their similarity to the HLA of interest 

	In the result files, all the p-values are ordered the same as the TCRs in the file TCR.txt. 

knn_split, 
	for those HLA with sample size between 40 and 70,
	it is slow to run

