library(pROC)
library(Matrix)

# this HLA_index is the HLA index in the dist_matrix
# For the HLA-I we didnt have this problem because all the rows in the matrix have its corresponding entry in the HLA matrix.
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
HLA_no_NA = HLA[,!is.na(HLA[216,])]
HLA_only = HLA_no_NA[1:215,]
CMV = HLA_no_NA[216,]
similar_amino = readRDS("../../data/correlation_matrices.rds")$aa



neigbor4_weights_mat = matrix(nrow = nrow(similar_amino), ncol=5)
neigbor4_index_mat = matrix(nrow = nrow(similar_amino), ncol=5)
threshold = 0.9




for (i in 1:nrow(similar_amino)) {
  NN_res = HLA_KNN(i, similar_amino, 5, HLA_no_NA, threshold)
  neigbor4_weights_mat[i,] = NN_res$neigbor_weights
  neigbor4_index_mat[i, ] = NN_res$neigbhor_index
  
}
####

output_dir = "../intermediate_files"
  
write.csv(neigbor4_weights_mat, 
          file.path(output_dir, "HLA_I_weights_mat.csv"), 
          row.names = F)

# index-1 for python format
write.csv(neigbor4_index_mat-1, 
          file.path(output_dir, "HLA_I_neighbor_index_mat.csv"), 
          row.names = F)

sessionInfo()
q(save = "no")
