

HLA_KNN = function(HLA_index, dist_mat, k, HLA, similar_threshold)
{
  test_weights = c(1, sort(dist_mat[HLA_index,], decreasing = T)[2:k])
  index_to_drop = which(test_weights < similar_threshold)
  
  test_weights_index = order(dist_mat[HLA_index,],decreasing = T )[1:k]
  HLA_at = which(test_weights_index ==  HLA_index)
  if (HLA_at != 1) {
    test_weights_index = test_weights_index[-HLA_at]
    test_weights_index = c(HLA_index, test_weights_index)
  }
  
  neigbhor_names = rownames(dist_mat)[test_weights_index]
  neigbhor_index_HLA_mat = match(neigbhor_names, row.names(HLA))
  test_weights[index_to_drop] = NA
  neigbhor_index_HLA_mat[index_to_drop] =NA
  
  return(list(neigbor_weights = test_weights, neigbhor_index = neigbhor_index_HLA_mat  ) )
}



HLA  = readRDS("../../data/DeWitt_2018/HLA_v2_CMV_data.rds")
HLA_no_NA = HLA[ ,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
CMV = HLA_no_NA[216,]
dat = readRDS("../../data/HLA_II_similarity_matrix_aa.rds")


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

# this is python stype 0-indexed HLA_II indexes
HLA_II_index=read.csv("../intermediate_files/HLA_II_index.csv",header = F)
dim(HLA_II_index)
head(HLA_II_index)

# this is R stype 1-indexed HLA_II indexes
HLA_II_index_verify = sort((   match(rownames(dat) , rownames(HLA_only))  )[!is.na(match(rownames(dat) , rownames(HLA_only)))])
HLA_II_index_verify = c(HLA_II_index_verify,c(216,217,218,219,220))
rownames(split_HLA)[HLA_II_index_verify]

# it verified that indexes from these two approached exactly match
table(HLA_II_index$V1+1 == HLA_II_index_verify, useNA="ifany")

write.csv(rownames(split_HLA)[HLA_II_index_verify],
          "../intermediate_files/HLA_II_names.csv",
          row.names = F)



dat_index = match(rownames(split_HLA)[HLA_II_index_verify], rownames(dat))
new_dat = dat[dat_index, dat_index]

neigbor4_weights_mat = matrix(nrow = length(dat_index), ncol=5)
neigbor4_index_mat = matrix(nrow = length(dat_index), ncol=5)
threshold = 0.9
  
for (i in 1:135) {
  print(i)
  NN_res = HLA_KNN(i, new_dat, 5, split_HLA, threshold)
  neigbor4_weights_mat[i,] = NN_res$neigbor_weights
  neigbor4_index_mat[i, ] = NN_res$neigbhor_index
}

## save tables out

output_dir = "../intermediate_files"

write.csv(neigbor4_weights_mat, 
          file.path(output_dir, "HLA_II_weights_mat.csv"),
          row.names = F)

# index-1 for python format
write.csv(neigbor4_index_mat-1, 
          file.path(output_dir, "HLA_II_neighbor_index_mat.csv"),
          row.names = F)


sessionInfo()
q(save = "no")




