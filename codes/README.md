
## HLA agnostic model

Run HLA agnostic model for each HLA from one training/test split:

	HLA_agnostic_model

For the HLAs with frequency in [40, 70], get predictions based on HLA agnostic models run on 10 training/test splits

	HLA_agnostic_model_split

## HLA specific model

For HLA-I alleles, HLA-II alleles, and the HLAs with frequency in [40, 70]. 

Also part of code files related to permutation tests and pvalues between TCR and HLA related to revision

	HLA_specific_model_HLA_I
	HLA_specific_model_HLA_II
	HLA_specific_model_split

## combined model

combined HLA model, for all HLAs, and the HLAs with frequency in [40, 70]

	combined_HLA_model
	combined_HLA_model_split

## KNN model

for HLA-I alleles, HLA-II alleles, and the HLAs with frequency in [40, 70]

	weighted_HLA_I
	weighted_HLA_II
	weighted_HLA_split

## HLA CMV association

Association between HLA and CMV

	HLA_CMV_association

## Code for data

Code files to generate the intermediate data files for running models

	code_for_data

## Code for tables

Code files to generate tables

	code_for_tables

## Code for revision

Code files related to certain discussions in revision

	code_for_revision

	





