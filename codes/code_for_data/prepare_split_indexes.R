library(dplyr)
set.seed(611)
HLA  = readRDS("../../data/DeWitt_2018/HLA_v2_CMV_data.rds")
HLA_no_NA = HLA[,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
CMV = HLA_no_NA[216,]

train_index_mat = matrix(nrow=10, ncol=512)
test_index_mat = matrix(nrow = 10, ncol = 129)
index_to_sample = 1:length(CMV)

for (i in 1:10) {
  train_index_mat[i,] =  sort(sample(index_to_sample, 512))
  test_index_mat[i,] = setdiff(index_to_sample, train_index_mat[i,])
}
#train_index_mat[1,]

HLA_I_index = read.csv("../../codes/intermediate_files/HLA_I_index.csv", header = F)$V1+1
HLA_II_index = read.csv("../../codes/intermediate_files/HLA_II_index.csv", header = F)$V1+1

# train must have 32, test must have 8
# all the HLAs in the HLA_matrix is either HLA_I or HLA_II
dim(HLA_no_NA)
subj_counts = rowSums(HLA_no_NA[-216,], na.rm = T)
#which(subj_counts >= 40 & subj_counts <= 70)
mid_HLAs =  subj_counts[subj_counts >= 40 & subj_counts <= 70]
mid_HLA_indexes = match(names(mid_HLAs), rownames(HLA_no_NA))

length(subj_counts[subj_counts < 40])
length(subj_counts[subj_counts < 30])

mid_HLA_indexes

# what I need now is a function that sort of filters it

HLA_mid_or_not_split = function(train_index, test_index, HLA_index, HLA_mat)
{
  train_sum = sum(HLA_mat[HLA_index, train_index], na.rm = T)
  test_sum = sum(HLA_mat[HLA_index, test_index], na.rm = T)
  
  if (train_sum >=32 & test_sum >= 8) {
    return(1)
  }
  return(0)
}

res_mat = matrix(nrow= length(mid_HLA_indexes), ncol = 10)


for (index in 1:length(mid_HLA_indexes)) {
  res_vec = numeric(10)
  for (i in 1:10) {
    res_vec[i] = HLA_mid_or_not_split(train_index_mat[i,], test_index_mat[i,], mid_HLA_indexes[index] , HLA_no_NA)
  }
  res_mat[index,] = res_vec
}

res_mat

mid_HLA_indexes

# account for haplotype HLA-DRDQ*13:01_01:03_06:03
# among the splits for it
# HLA-DRB*13:01 takes row 89 (R index)
# in the python code split_knn.py, HLA-DQAB*01:03_06:03 takes row 219 (R index)
mid_HLA_indexes = c(mid_HLA_indexes, 219)
mid_HLA_indexes
length(mid_HLA_indexes)

# since they form a haplotype, they appear exactly in the same set of individuals
rownames(HLA_no_NA)[89]
which(mid_HLA_indexes==89)
res_mat = rbind(res_mat, res_mat[which(mid_HLA_indexes==89), ])

write.table(res_mat, "../../codes/intermediate_files/split_mat.csv", row.names = F, sep=",",col.names = F)
write.table(mid_HLA_indexes-1, "../../codes/intermediate_files/mid_HLA_index.csv", row.names = F, col.names = F)
write.table(train_index_mat-1, "../../codes/intermediate_files/split_train_index.csv", row.names = F, sep=",",col.names = F)
write.table(test_index_mat-1, "../../codes/intermediate_files/split_test_index.csv", row.names = F, sep=",",col.names = F)


sessionInfo()
q(save = "no")

