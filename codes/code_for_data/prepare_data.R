library(Matrix)
library(pROC)

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



