library(dplyr)

HLA  = readRDS("../conditional_TCR_prediction/data/DeWitt_2018/HLA_v2_CMV_data.rds")
HLA_no_NA = HLA[ ,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
CMV = HLA_no_NA[216,]
dat = readRDS("../conditional_TCR_prediction/data/HLA_II_similarity_matrix_aa.rds")

#These are the HLAs that has been lumped together, now we are breaking them.

exceptions = c("HLA-DRDQ*03:01_05:01_02:01","HLA-DRDQ*09:01_03:02_03:03",
               "HLA-DRDQ*10:01_01:05_05:01","HLA-DRDQ*13:01_01:03_06:03", "HLA-DRDQ*15:01_01:02_06:02")
exception_DRB1_names = c("HLA-DRB1*03:01", "HLA-DRB1*09:01", "HLA-DRB1*10:01",
                         "HLA-DRB1*13:01", "HLA-DRB1*15:01")
exception_DQAB_names = c("HLA-DQAB*05:01_02:01", "HLA-DQAB*03:02_03:03", "HLA-DQAB*01:05_05:01",
                         "HLA-DQAB*01:03_06:03", "HLA-DQAB*01:02_06:02")
exception_index = match(exceptions, rownames(HLA_only))
exception_index

newHLA_names = rownames(HLA_only)
newdata = HLA_only[exception_index,]
rownames(newdata) = exception_DQAB_names
newHLA_names[exception_index] = exception_DRB1_names
rownames(HLA_only) = newHLA_names
split_HLA = rbind(HLA_only,newdata)

# This extracts all the HLAs that are both in the correlation matrix as well as the HLA matrix.
HLA_II_index = sort((match(rownames(dat) , rownames(HLA_only)))[!is.na(match(rownames(dat) , rownames(HLA_only)))])

# This concates the extra HLAs
HLA_II_index = c(HLA_II_index,c(216,217,218,219,220))
rownames(split_HLA)[HLA_II_index]

# Outupts the name of the HLA-IIs
write.csv(rownames(split_HLA)[HLA_II_index],"HLA_II_names.csv",row.names = F)


# This outputs the index for python, subtract the index by 1 since python index starts at 0.
write.table(HLA_II_index-1,"HLA_II_index.csv",row.names = F, col.names = FALSE)

# These are the HLA-II alleles that have higher score in specific model
hla_ii_higher = c("HLA-DRB1*07:01", "HLA-DRB1*03:01", 
                  "HLA-DQAB*01:01_05:01", "HLA-DQAB*02:01_02:02", 
                  "HLA-DQAB*05:01_02:01")

# get the 0-indexed relative indices for these alleles
which(rownames(split_HLA)[HLA_II_index]%in%hla_ii_higher) - 1

rownames(split_HLA)[HLA_II_index][which(rownames(split_HLA)[HLA_II_index]%in%hla_ii_higher)]
