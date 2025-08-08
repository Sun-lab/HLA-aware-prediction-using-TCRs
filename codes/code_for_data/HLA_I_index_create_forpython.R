library(dplyr)

HLA  = readRDS("../../data/DeWitt_2018/HLA_v2_CMV_data.rds")
HLA_no_NA = HLA[ ,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
CMV = HLA_no_NA[216,]
dat = readRDS("../../data/correlation_matrices.rds")$aa

match(rownames(dat) , rownames(HLA_only))
HLA_I_index = sort((match(rownames(dat) , rownames(HLA_only)))[!is.na(match(rownames(dat) , rownames(HLA_only)))])

# names
write.csv(rownames(HLA_only)[HLA_I_index],"../intermediate_files/HLA_I_names.csv",row.names = F)

# This outputs the index for python, subtract the index by 1 since python index starts at 0.
write.table(HLA_I_index-1,"../intermediate_files/HLA_I_index.csv",row.names = F, col.names=FALSE)


sessionInfo()
q(save = "no")