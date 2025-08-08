library(dplyr)

HLA  = readRDS("../conditional_TCR_prediction/data/DeWitt_2018/HLA_v2_CMV_data.rds")
HLA_no_NA = HLA[,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
CMV = HLA_no_NA[216,]

# all the HLAs in the HLA_matrix is either HLA_I or HLA_II
# HLA_no_NA
subj_counts = rowSums(HLA_no_NA[-216,], na.rm = T)
#which(subj_counts >= 40 & subj_counts <= 70)
mid_HLAs =  subj_counts[subj_counts >= 40 & subj_counts <= 70]
mid_HLA_indexes = match(names(mid_HLAs), rownames(HLA_no_NA))

length(subj_counts[subj_counts < 40])
length(subj_counts[subj_counts < 30])

# according to the content of mid_HLA_indexes, the HLA with python index 88 (R index 89)
# is among those with frequency between [40, 70]
# and it is also one of the haplotypes "HLA-DRDQ*13:01_01:03_06:03"
# so HLA-DQAB01:03_06:03 also has frequency in [40 70]
# and should have the same AUC as HLA-DRB*13:01
# from Table S1, it seems not the case
mid_HLA_indexes


# direct load the saved out full HLA matrix with split haplotypes and NAs kept
HLA_only = read.csv(file.path("../intermediate_files", "complete_HLA_NAs_kept.csv"), 
                    header=TRUE)
dim(HLA_only)

subj_counts = rowSums(HLA_only, na.rm = T)
high_HLA_indexes = which(subj_counts > 70)
mid_HLA_indexes = which((subj_counts >= 40) & (subj_counts <= 70))

write.table(high_HLA_indexes-1, 
            "../intermediate_files/high_HLA_index_complete.csv", 
            sep = ",", row.names = F, col.names = F)

write.table(mid_HLA_indexes-1, 
            "../intermediate_files/mid_HLA_index.csv", 
            sep = ",", row.names = F, col.names = F)


sessionInfo()
q(save = "no")