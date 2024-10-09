library(Matrix)
library(pROC)

setwd("./")
TCR_path = "../../data/Emerson_2017/TCR_data.rds"
HLA_path = "../../data/DeWitt_2018/HLA_v2_CMV_data.rds"
dest_TCR_path = "../intermediate_files/TCR.txt"
dest_HLA_path = "../intermediate_files/HLA.txt"
dest_CMV_path = "../intermediate_files/CMV.txt"


# HLA & CMV 
HLA  = readRDS(HLA_path)
HLA_no_NA = HLA[,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
CMV = HLA_no_NA[216,]

# missing info treated as not having
HLA_only[is.na(HLA_only)] <- 0

# TCR
TCR = readRDS(TCR_path)
TCR_no_NA = TCR[,!is.na(HLA[216,])]

HLA_sparse = Matrix(HLA_only, sparse = TRUE)
writeMM(TCR_no_NA, file = dest_TCR_path)
writeMM(HLA_sparse, file = dest_HLA_path)
write.csv(CMV, file = dest_CMV_path)

# Complete HLA that contains the couple HLA-IIs

dat = readRDS("../../data/HLA_II_similarity_matrix_aa.rds")
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
saveRDS(split_HLA, "../intermediate_files/complete_HLA.rds")


