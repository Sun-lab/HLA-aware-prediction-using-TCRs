# in this file, the HLA matrix is processed to keep the NA encoding
# to prepare for running association test between HLA and TCRs

library(Matrix)

# the current HLA processing procedure is to only drop the individuals that doesn't have CMV status
# need to keep the NA encoding for individuals who do not have HLA typing results for the location

HLA  = readRDS("../conditional_TCR_prediction/data/DeWitt_2018/HLA_v2_CMV_data.rds")
HLA_no_NA = HLA[,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
# need to keep NAs for the missing info 

# write out an HLA file for the purpose of association tests either with TCR or CMV
dat = readRDS("../conditional_TCR_prediction/data/HLA_II_similarity_matrix_aa.rds")
exceptions = c("HLA-DRDQ*03:01_05:01_02:01","HLA-DRDQ*09:01_03:02_03:03",
               "HLA-DRDQ*10:01_01:05_05:01","HLA-DRDQ*13:01_01:03_06:03", "HLA-DRDQ*15:01_01:02_06:02")

exception_DRB1_names = c("HLA-DRB1*03:01", "HLA-DRB1*09:01", "HLA-DRB1*10:01",
                         "HLA-DRB1*13:01", "HLA-DRB1*15:01")

exception_DQAB_names = c("HLA-DQAB*05:01_02:01", "HLA-DQAB*03:02_03:03", "HLA-DQAB*01:05_05:01",
                         "HLA-DQAB*01:03_06:03", "HLA-DQAB*01:02_06:02")

exception_index = match(exceptions, rownames(HLA_only))
newHLA_names = rownames(HLA_only)
newdata = HLA_only[exception_index,]
rownames(newdata) = exception_DQAB_names
newHLA_names[exception_index] = exception_DRB1_names
rownames(HLA_only) = newHLA_names
split_HLA = rbind(HLA_only,newdata)
saveRDS(split_HLA, "../intermediate_files/complete_HLA_NAs_kept.rds")

df_split_HLA = as.data.frame(split_HLA)

df_split_HLA_rownames = data.frame(HLA_name = row.names(df_split_HLA))

write.csv(df_split_HLA_rownames, 
          "../intermediate_files/complete_HLA_rownames.csv", 
          row.names = FALSE)

write.csv(df_split_HLA, 
          "../intermediate_files/complete_HLA_NAs_kept.csv", 
          row.names = FALSE)

sessionInfo()
q(save = "no")
