# generate Table S1

# three sheets
# first sheet: explanation for each column name
# second sheet: HLA-I alleles results
# third sheet: HLA-II alleles results

# two versions of tables:
# one with only the HLAs having frequency at least 40
# the other one with all HLAs

# for HLAs with frequency above 70 or below 40, 
#    the values are from only one split
# for HLAs with frequency in [40, 70], 
#    the values are averaged from the valid splits out of the 10 splits

# columns to have:
# 1. HLA
# 2. HLA_agnostic_AUC
# 3. HLA_specific_AUC
# 4. combined_AUC
# 5. HLA_specific_training_size
# 6. Number of TCRs significant under HLA-agnostic model
# 7. Number of TCRs significant under HLA-specific model
# 8. Size of intersection between HLA-specific significant TCRs and HLA_agnostic significant TCRs
# 9. Number of TCRs that are significant under HLA-agnostic model, but not under HLA-specific model
# 10. Number of TCRs that are significant under HLA-specific model, but not under HLA-agnostic model
# 11. Proportion of intersection out of HLA_agnostic significant TCRs
# 12. Proportion of intersection out of HLA-specific significant TCRs
# 13. KNN_AUC
# 14. KNN_training_size
# 15. Number of TCRs significant under KNN model
# 16. Size of intersection between KNN model significant TCRs and HLA_agnostic significant TCRs
# 17. Number of TCRs that are significant under HLA-agnostic model, but not under KNN model
# 18. Number of TCRs that are significant under KNN model, but not under HLA-agnostic model
# 19. Proportion of intersection out of KNN significant TCRs
# 20. Proportion of intersection out of HLA_agnostic significant TCRs
# 21. HLA_group
# 22. HLA_name
# 23. HLA_frequency
# 24. run_10splits

library(openxlsx)


file_dir = "../intermediate_files"

pval_cutoff = 0.001

df_HLA = readRDS(file.path(file_dir, "complete_HLA_NAs_kept.rds"))
dim(df_HLA)


# 1. HLA ---------------------------

col_HLA = rownames(df_HLA)
  
# the indexes are 0-indexed, following the format in Python
# mid_HLA_index are for the HLAs that have frequency between [40, 70]
#  these HLAs go through 10 splits
df_HLA_I_index = read.csv(file.path(file_dir, "HLA_I_index.csv"), header=FALSE)
df_HLA_II_index = read.csv(file.path(file_dir, "HLA_II_index.csv"), header=FALSE)
df_HLA_mid = read.csv(file.path(file_dir, "mid_HLA_index.csv"), header=FALSE)

HLA_I_index = df_HLA_I_index$V1
HLA_II_index = df_HLA_II_index$V1
HLA_mid = df_HLA_mid$V1

# 2. HLA_agnostic_AUC ---------------------------

agnostic_dir = "../HLA_agnostic_model/results"
agnostic_split_dir = "../HLA_agnostic_model_split/results"

agnostic_auc = read.csv(file.path(agnostic_dir, "specific_aucs.csv"), 
                        header=TRUE)
agnostic_auc = agnostic_auc$AUC

stopifnot(sum(agnostic_auc<0, na.rm=TRUE)==0)

agnostic_split_auc_matrix = NULL

for (split_i in 0:9){
  df_agnostic_split_i = read.csv(file.path(agnostic_split_dir, 
                                           paste0("aucs_split", as.character(split_i), ".csv")), 
                                 header=TRUE)
  stopifnot(sum(df_agnostic_split_i$AUC<0, na.rm=TRUE)==0)
  agnostic_split_auc_matrix = cbind(agnostic_split_auc_matrix, df_agnostic_split_i$AUC)
}

stopifnot(all(HLA_mid==df_agnostic_split_i$HLA_index))

agnostic_split_auc_mean = rowMeans(agnostic_split_auc_matrix, na.rm=TRUE)

col_HLA_agnostic_AUC = rep(NA, nrow(df_HLA))

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    col_HLA_agnostic_AUC[i+1] = agnostic_split_auc_mean[which_row]
  }else{
    col_HLA_agnostic_AUC[i+1] = agnostic_auc[i+1]
  }
}

# 2. HLA_specific_AUC ---------------------------

specific_split_dir = "../HLA_specific_model_split/results"
  
specific_split_auc_matrix = NULL

for (hla_i in HLA_mid){
  df_specific_hla_i = read.csv(file.path(specific_split_dir, 
                                           paste0("res_aucs_hla_index_", as.character(hla_i)),
                                           paste0("aucs_hla_index_", as.character(hla_i), ".csv")), 
                                 header=TRUE)
  stopifnot(sum(df_specific_hla_i$AUC<0, na.rm=TRUE)==0)
  specific_split_auc_matrix = rbind(specific_split_auc_matrix, df_specific_hla_i$AUC)
}

specific_split_auc_mean = rowMeans(specific_split_auc_matrix, na.rm=TRUE)


col_HLA_specific_AUC = rep(NA, nrow(df_HLA))

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    col_HLA_specific_AUC[i+1] = specific_split_auc_mean[which_row]
  }else{
    if (i%in%HLA_I_index){
      specific_dir = "../HLA_specific_model_HLA_I/results"
    }else{
      stopifnot(i%in%HLA_II_index)
      specific_dir = "../HLA_specific_model_HLA_II/results"      
    }
    df_specific_i = read.csv(file.path(specific_dir, "aucs", 
                                       paste0("aucs_", as.character(i), ".csv")), 
                             header=TRUE)  
    stopifnot(sum(df_specific_i$AUC<0, na.rm=TRUE)==0)
    auc_i = df_specific_i$AUC[2] # the location corresponding to pvalue cutoff 0.001
    col_HLA_specific_AUC[i+1] = auc_i
  }
}

# 4. combined_AUC ---------------------------

combined_dir = "../combined_HLA_model/results"
combined_split_dir = "../combined_HLA_model_split/results"

combined_auc = read.csv(file.path(combined_dir, "combined_aucs.csv"), 
                        header=TRUE)
combined_auc = combined_auc$AUC

stopifnot(sum(combined_auc<0, na.rm=TRUE)==0)

combined_split_auc_matrix = NULL

for (hla_i in HLA_mid){
  df_combined_hla_i = read.csv(file.path(combined_split_dir, 
                                         paste0("combined_split_aucs_hla_index_", as.character(hla_i), ".csv")), 
                               header=TRUE)
  stopifnot(sum(df_combined_hla_i$AUC<0, na.rm=TRUE)==0)
  combined_split_auc_matrix = rbind(combined_split_auc_matrix, df_combined_hla_i$AUC)
}

combined_split_auc_mean = rowMeans(combined_split_auc_matrix, na.rm=TRUE)


col_combined_AUC = rep(NA, nrow(df_HLA))

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    col_combined_AUC[i+1] = combined_split_auc_mean[which_row]
  }else{
    col_combined_AUC[i+1] = combined_auc[i+1]
  }
}

# 5. HLA_specific_training_size ---------------------------
# 23. HLA_frequency
# for HLAs with frequency between [40, 70]
# take the average of the training sizes from all valid splits

df_train_index = read.csv(file.path(file_dir, "train_index.csv"), header=FALSE)
df_split_mat = read.csv(file.path(file_dir, "split_mat.csv"), header=FALSE)
df_split_train_index = read.csv(file.path(file_dir, "split_train_index.csv"), header=FALSE)

train_index = df_train_index$V1
split_mat = as.matrix(df_split_mat)
split_train_index = as.matrix(df_split_train_index)

col_HLA_specific_training_size = rep(NA, nrow(df_HLA))

HLA_mat = as.matrix(df_HLA)

HLA_mat_train = HLA_mat[, train_index+1]
HLA_training_sizes = rowSums(HLA_mat_train, na.rm=TRUE) # from major split

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    # only use the splits that are valid under the given HLA
    valid_splits = which(as.numeric(split_mat[which_row,])==1)
    train_size_splits = NULL
    for (split_i in valid_splits){
      split_train_index_i = as.numeric(split_train_index[split_i, ])
      HLA_mat_train_split_i = HLA_mat[, split_train_index_i+1]
      HLA_training_sizes_split_i = as.numeric(rowSums(HLA_mat_train_split_i, na.rm=TRUE))
      train_size_splits = c(train_size_splits, HLA_training_sizes_split_i[i+1])
    }
    col_HLA_specific_training_size[i+1] = mean(train_size_splits)
  }else{
    col_HLA_specific_training_size[i+1] = HLA_training_sizes[i+1]
  }
}

# 6. Number of TCRs significant under HLA-agnostic model ---------------------------
# 7. Number of TCRs significant under HLA-specific model ---------------------------
# 8. Size of intersection between HLA-specific significant TCRs and HLA_agnostic significant TCRs
# 9. Number of TCRs that are significant under HLA-agnostic model, but not under HLA-specific model
# 10. Number of TCRs that are significant under HLA-specific model, but not under HLA-agnostic model
# 11. Proportion of intersection out of HLA_agnostic significant TCRs
# 12. Proportion of intersection out of HLA-specific significant TCRs

# prepare TCRs significant under HLA-agnostic model

agnostic_dir = "../HLA_agnostic_model/results"
agnostic_split_dir = "../HLA_agnostic_model_split/results"

agnostic_pvals = read.csv(file.path(agnostic_dir, "pvalues.csv"), 
                          header=TRUE)

tcr_index_agnostic = which(agnostic_pvals$pval<=pval_cutoff) # for sign TCR set operations

#n_sign_agnostic = sum(agnostic_pvals$pval<=pval_cutoff, na.rm=TRUE)

tcr_index_agnostic_split_list = list() # for sign TCR set operations

for (split_i in 0:9){
  agnostic_pvals_i = read.csv(file.path(agnostic_split_dir, 
                                        paste0("pvalues_split", as.character(split_i), ".csv")), 
                              header=TRUE)
  tcr_index_agnostic_split_list[[paste0("split", as.character(split_i))]] = 
    which(agnostic_pvals_i$pval<=pval_cutoff)

}
  
# prepare TCRs significant under HLA-specific model

tcr_index_specific_list = list() # for sign TCR set operations
tcr_index_specific_split_list = list() # for sign TCR set operations

specific_split_dir = "../HLA_specific_model_split/results"

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    # only use the splits that are valid under the given HLA
    valid_splits = which(as.numeric(split_mat[which_row,])==1)
    cur_sign_tcr_list = list()
    for (split_i in valid_splits){
      specific_pvals_split_i = read.csv(file.path(specific_split_dir, 
                                                  paste0("pvals_hla_index_", as.character(i)), 
                                                  paste0("pvalues_hla_index_", as.character(i), 
                                                         "_", as.character(split_i-1),".csv")), 
                                        header=TRUE)
      cur_sign_tcr_list[[sprintf("split%d", (split_i-1))]] = 
        which(specific_pvals_split_i$pval<=pval_cutoff)
    }
    tcr_index_specific_split_list[[sprintf("hla_%d", i)]] = cur_sign_tcr_list
  }else{
    if (i%in%HLA_I_index){
      specific_dir = "../HLA_specific_model_HLA_I/results"
    }else{
      stopifnot(i%in%HLA_II_index)
      specific_dir = "../HLA_specific_model_HLA_II/results"      
    }
    specific_pvals_i = read.csv(file.path(specific_dir, "pvals", 
                                          sprintf("pvalues_hla_index_%d.csv", i)), 
                                header=TRUE)  
    tcr_index_specific_list[[sprintf("hla_%d", i)]] = 
      which(specific_pvals_i$pval<=pval_cutoff)
  }
}

length(tcr_index_specific_list)
length(tcr_index_specific_split_list)

# get values for columns
col_n_sign_agnostic = rep(NA, nrow(df_HLA))
col_n_sign_specific = rep(NA, nrow(df_HLA))
col_n_intersect_agnostic_specific = rep(NA, nrow(df_HLA))
col_n_agnostic_only = rep(NA, nrow(df_HLA))
col_n_specific_only = rep(NA, nrow(df_HLA))
col_prop_divided_by_agnostic = rep(NA, nrow(df_HLA))
col_prop_divided_by_specific = rep(NA, nrow(df_HLA))

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    # only use the splits that are valid under the given HLA
    valid_splits = which(as.numeric(split_mat[which_row,])==1)
    cur_n_sign_agnostic_splits = NULL
    cur_n_sign_specific_splits = NULL
    cur_n_intersect_agnostic_specific_splits = NULL
    cur_n_agnostic_only_splits = NULL
    cur_n_specific_only_splits = NULL
    cur_prop_divided_by_agnostic_splits = NULL
    cur_prop_divided_by_specific_splits = NULL

    for (split_i in valid_splits){

      # sign agnostic from split_i
      cur_tcrs_agnostic_split = tcr_index_agnostic_split_list[[sprintf("split%d", (split_i-1))]]
      # sign specific from split_i
      cur_tcrs_specific_split = tcr_index_specific_split_list[[sprintf("hla_%d", i)]][[sprintf("split%d", (split_i-1))]]
      # intersect
      cur_intersect = intersect(cur_tcrs_agnostic_split, 
                                cur_tcrs_specific_split)
      
      cur_n_sign_agnostic_splits = c(cur_n_sign_agnostic_splits, 
                                     length(cur_tcrs_agnostic_split))
      cur_n_sign_specific_splits = c(cur_n_sign_specific_splits, 
                                     length(cur_tcrs_specific_split))
      cur_n_intersect_agnostic_specific_splits = c(cur_n_intersect_agnostic_specific_splits, 
                                                   length(cur_intersect))
      cur_n_agnostic_only_splits = c(cur_n_agnostic_only_splits, 
                                     length(setdiff(cur_tcrs_agnostic_split, cur_intersect)))
      cur_n_specific_only_splits = c(cur_n_specific_only_splits, 
                                     length(setdiff(cur_tcrs_specific_split, cur_intersect)))
      cur_prop_divided_by_agnostic_splits = c(cur_prop_divided_by_agnostic_splits, 
                                              length(cur_intersect)/length(cur_tcrs_agnostic_split))
      cur_prop_divided_by_specific_splits = c(cur_prop_divided_by_specific_splits, 
                                              length(cur_intersect)/length(cur_tcrs_specific_split))

    }
    
    col_n_sign_agnostic[i+1] = mean(cur_n_sign_agnostic_splits)
    col_n_sign_specific[i+1] = mean(cur_n_sign_specific_splits)
    col_n_intersect_agnostic_specific[i+1] = mean(cur_n_intersect_agnostic_specific_splits)
    col_n_agnostic_only[i+1] = mean(cur_n_agnostic_only_splits)
    col_n_specific_only[i+1] = mean(cur_n_specific_only_splits)
    col_prop_divided_by_agnostic[i+1] = mean(cur_prop_divided_by_agnostic_splits)
    col_prop_divided_by_specific[i+1] = mean(cur_prop_divided_by_specific_splits, na.rm=TRUE)
    
  }else{
 
    cur_tcrs_specific = tcr_index_specific_list[[sprintf("hla_%d", i)]]   
    
    cur_intersect = intersect(tcr_index_agnostic, 
                              cur_tcrs_specific)
    
    n_sign_agnostic = length(tcr_index_agnostic)
    
    col_n_sign_agnostic[i+1] = n_sign_agnostic
    col_n_sign_specific[i+1] = length(cur_tcrs_specific)
    col_n_intersect_agnostic_specific[i+1] = length(cur_intersect)
    col_n_agnostic_only[i+1] = length(setdiff(tcr_index_agnostic, cur_intersect))
    col_n_specific_only[i+1] = length(setdiff(cur_tcrs_specific, cur_intersect))
    col_prop_divided_by_agnostic[i+1] = length(cur_intersect)/n_sign_agnostic
    col_prop_divided_by_specific[i+1] = length(cur_intersect)/length(cur_tcrs_specific)

  }
}
    

# 13. KNN_AUC ---------------------------

weighted_split_dir = "../weighted_HLA_split/results"

col_HLA_weighted_AUC = rep(NA, nrow(df_HLA))

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    # only use the splits that are valid under the given HLA
    valid_splits = which(as.numeric(split_mat[which_row,])==1)   
    cur_auc_splits = NULL
    for (split_i in valid_splits){    
      df_auc_hla_split = read.csv(file.path(weighted_split_dir, 
                                            sprintf("res_aucs_%d", i),
                                            sprintf("auc_%d_%d.csv", i, (split_i-1))), 
                                  header=TRUE)
      stopifnot(nrow(df_auc_hla_split)==1)
      stopifnot(sum(df_auc_hla_split$AUC<0, na.rm=TRUE)==0)
      cur_auc_splits = c(cur_auc_splits, df_auc_hla_split$AUC)
    }
    col_HLA_weighted_AUC[i+1] = mean(cur_auc_splits, na.rm=TRUE)
  }else{
    if (i%in%HLA_I_index){
      weighted_dir = "../weighted_HLA_I/results"
    }else{
      stopifnot(i%in%HLA_II_index)
      weighted_dir = "../weighted_HLA_II/results"      
    }
    df_weighted_i = read.csv(file.path(weighted_dir, "res_aucs", 
                                       sprintf("auc_hla_index_%d.csv", i)), 
                             header=TRUE)  
    stopifnot(sum(df_weighted_i$AUC<0, na.rm=TRUE)==0)
    auc_i = df_weighted_i$AUC[1] # the location corresponding to the center HLA
    col_HLA_weighted_AUC[i+1] = auc_i
  }
}






# 14. KNN_training_size ---------------------------

col_HLA_weighted_training_size = rep(NA, nrow(df_HLA))

df_HLA_I_neighbor = read.csv("../intermediate_files/HLA_I_neighbor_index_mat.csv", 
                             header = TRUE)
df_HLA_II_neighbor = read.csv("../intermediate_files/HLA_II_neighbor_index_mat.csv", 
                             header = TRUE)

dim(df_HLA_I_neighbor)
dim(df_HLA_II_neighbor)
                   
HLA_I_neighbor = as.matrix(df_HLA_I_neighbor)
HLA_II_neighbor = as.matrix(df_HLA_II_neighbor)

# train_index = df_train_index$V1
# split_mat = as.matrix(df_split_mat)
# split_train_index = as.matrix(df_split_train_index)

# HLA_mat = as.matrix(df_HLA)

HLA_mat_train = HLA_mat[, train_index+1]


for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    # only use the splits that are valid under the given HLA
    valid_splits = which(as.numeric(split_mat[which_row,])==1)
    train_size_splits = NULL
    for (split_i in valid_splits){
      split_train_index_i = as.numeric(split_train_index[split_i, ])
      HLA_mat_train_split_i = HLA_mat[, split_train_index_i+1]
      # count the training individuals with at least one of the HLAs in the neighborhood
      if (i%in%HLA_I_index){
        which_i = which(HLA_I_index==i)
        neighbors_i = as.numeric(HLA_I_neighbor[which_i,])
      }else{
        stopifnot(i%in%HLA_II_index)
        which_i = which(HLA_II_index==i)
        neighbors_i = as.numeric(HLA_II_neighbor[which_i,]) 
      }
      neighbors_i = neighbors_i[!is.na(neighbors_i)] 
      
      HLA_mat_slice = HLA_mat_train_split_i[(neighbors_i+1), , drop = FALSE]
      
      train_size_splits = c(train_size_splits, 
                            sum(colSums(HLA_mat_slice, na.rm=TRUE)>0))
      
    }
    col_HLA_weighted_training_size[i+1] = mean(train_size_splits)
  }else{
    if (i%in%HLA_I_index){
      which_i = which(HLA_I_index==i)
      neighbors_i = as.numeric(HLA_I_neighbor[which_i,])
    }else{
      stopifnot(i%in%HLA_II_index)
      which_i = which(HLA_II_index==i)
      neighbors_i = as.numeric(HLA_II_neighbor[which_i,]) 
    }
    neighbors_i = neighbors_i[!is.na(neighbors_i)] 
    
    HLA_mat_slice = HLA_mat_train[(neighbors_i+1), , drop = FALSE]    
    
    col_HLA_weighted_training_size[i+1] = 
      sum(colSums(HLA_mat_slice, na.rm=TRUE)>0)
  }
}



# 15. Number of TCRs significant under KNN model ---------------------------
# 16. Size of intersection between KNN model significant TCRs and HLA_agnostic significant TCRs
# 17. Number of TCRs that are significant under HLA-agnostic model, but not under KNN model
# 18. Number of TCRs that are significant under KNN model, but not under HLA-agnostic model
# 19. Proportion of intersection out of KNN significant TCRs
# 20. Proportion of intersection out of HLA_agnostic significant TCRs


# prepare TCRs significant under HLA-weighted model

tcr_index_weighted_list = list() # for sign TCR set operations
tcr_index_weighted_split_list = list() # for sign TCR set operations

weighted_split_dir = "../weighted_HLA_split/results"

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    # only use the splits that are valid under the given HLA
    valid_splits = which(as.numeric(split_mat[which_row,])==1)
    cur_sign_tcr_list = list()
    for (split_i in valid_splits){
      weighted_tcrs_split_i = read.csv(file.path(weighted_split_dir, 
                                                  sprintf("sig_indexes_%d", i), 
                                                  sprintf("sig_indexes_%d_%d.csv", i, (split_i-1))), 
                                        header=TRUE)
      cur_sign_tcr_list[[sprintf("split%d", (split_i-1))]] = 
        (weighted_tcrs_split_i$significant_TCR_index+1)
    }
    tcr_index_weighted_split_list[[sprintf("hla_%d", i)]] = cur_sign_tcr_list
  }else{
    if (i%in%HLA_I_index){
      weighted_dir = "../weighted_HLA_I/results"
    }else{
      stopifnot(i%in%HLA_II_index)
      weighted_dir = "../weighted_HLA_II/results"      
    }
    weighted_pvals_i = read.csv(file.path(weighted_dir, "pvalues", 
                                          sprintf("pvalues_hla_index_%d.csv", i)), 
                                header=TRUE)  
    tcr_index_weighted_list[[sprintf("hla_%d", i)]] = 
      which(weighted_pvals_i$pval<=pval_cutoff)
  }
}

length(tcr_index_weighted_list)
length(tcr_index_weighted_split_list)


# get values for columns
col_n_sign_weighted = rep(NA, nrow(df_HLA))
col_n_intersect_agnostic_weighted = rep(NA, nrow(df_HLA))
col_n_agnostic_not_in_weighted = rep(NA, nrow(df_HLA))
col_n_weighted_not_in_agnostic = rep(NA, nrow(df_HLA))
col_both_agnostic_weighted_divided_by_agnostic = rep(NA, nrow(df_HLA))
col_both_agnostic_weighted_divided_by_weighted = rep(NA, nrow(df_HLA))

for (i in 0:219){
  if (i%in%HLA_mid){
    which_row = which(HLA_mid==i)
    # only use the splits that are valid under the given HLA
    valid_splits = which(as.numeric(split_mat[which_row,])==1)

    cur_n_sign_weighted_splits = NULL
    cur_n_intersect_agnostic_weighted_splits = NULL
    cur_n_agnostic_only_splits = NULL
    cur_n_weighted_only_splits = NULL
    cur_prop_divided_by_agnostic_splits = NULL
    cur_prop_divided_by_weighted_splits = NULL
    
    for (split_i in valid_splits){
      
      # sign agnostic from split_i
      cur_tcrs_agnostic_split = tcr_index_agnostic_split_list[[sprintf("split%d", (split_i-1))]]
      # sign weighted from split_i
      cur_tcrs_weighted_split = tcr_index_weighted_split_list[[sprintf("hla_%d", i)]][[sprintf("split%d", (split_i-1))]]
      # intersect
      cur_intersect = intersect(cur_tcrs_agnostic_split, 
                                cur_tcrs_weighted_split)
      
      cur_n_sign_weighted_splits = c(cur_n_sign_weighted_splits, 
                                     length(cur_tcrs_weighted_split))
      cur_n_intersect_agnostic_weighted_splits = c(cur_n_intersect_agnostic_weighted_splits, 
                                                   length(cur_intersect))
      cur_n_agnostic_only_splits = c(cur_n_agnostic_only_splits, 
                                     length(setdiff(cur_tcrs_agnostic_split, cur_intersect)))
      cur_n_weighted_only_splits = c(cur_n_weighted_only_splits, 
                                     length(setdiff(cur_tcrs_weighted_split, cur_intersect)))
      cur_prop_divided_by_agnostic_splits = c(cur_prop_divided_by_agnostic_splits, 
                                              length(cur_intersect)/length(cur_tcrs_agnostic_split))
      cur_prop_divided_by_weighted_splits = c(cur_prop_divided_by_weighted_splits, 
                                              length(cur_intersect)/length(cur_tcrs_weighted_split))
      
    }
    
    col_n_sign_weighted[i+1] = mean(cur_n_sign_weighted_splits)
    col_n_intersect_agnostic_weighted[i+1] = mean(cur_n_intersect_agnostic_weighted_splits)
    col_n_agnostic_not_in_weighted[i+1] = mean(cur_n_agnostic_only_splits)
    col_n_weighted_not_in_agnostic[i+1] = mean(cur_n_weighted_only_splits)
    col_both_agnostic_weighted_divided_by_agnostic[i+1] = mean(cur_prop_divided_by_agnostic_splits)
    col_both_agnostic_weighted_divided_by_weighted[i+1] = mean(cur_prop_divided_by_weighted_splits, na.rm=TRUE)
    
  }else{
    
    cur_tcrs_weighted = tcr_index_weighted_list[[sprintf("hla_%d", i)]]   
    
    cur_intersect = intersect(tcr_index_agnostic, 
                              cur_tcrs_weighted)
    
    n_sign_agnostic = length(tcr_index_agnostic)
    
    col_n_sign_weighted[i+1] = length(cur_tcrs_weighted)
    col_n_intersect_agnostic_weighted[i+1] = length(cur_intersect)
    col_n_agnostic_not_in_weighted[i+1] = length(setdiff(tcr_index_agnostic, cur_intersect))
    col_n_weighted_not_in_agnostic[i+1] = length(setdiff(cur_tcrs_weighted, cur_intersect))
    col_both_agnostic_weighted_divided_by_agnostic[i+1] = length(cur_intersect)/n_sign_agnostic
    col_both_agnostic_weighted_divided_by_weighted[i+1] = length(cur_intersect)/length(cur_tcrs_weighted)
    
  }
}


# 21. HLA_group ---------------------------
# 22. HLA_name ---------------------------

col_HLA_name = substr(rownames(df_HLA), 5, nchar(rownames(df_HLA)))

HLA_types = sub("\\*.*", "", col_HLA_name)

col_HLA_group = substr(HLA_types, 1, 2)


# 23. HLA_frequency ---------------------------
# the frequency of HLAs among the 641 individuals with known CMV status (not being NA)
col_HLA_frequency = as.numeric(rowSums(HLA_mat, na.rm=TRUE))

# 24. run_10splits
col_run_10splits = rep(FALSE, nrow(df_HLA))
col_run_10splits[which((col_HLA_frequency>=40)&(col_HLA_frequency<70))] = TRUE

# 25. CMV_asso_pvalue ---------------------------
df_cmv_pval1 = read.csv("../HLA_CMV_association/results/HLA_I_CMV_asso_res.csv", 
                     header=TRUE)
df_cmv_pval2 = read.csv("../HLA_CMV_association/results/HLA_II_CMV_asso_res.csv", 
                     header=TRUE)
dim(df_cmv_pval1)
dim(df_cmv_pval2)

colnames(df_cmv_pval1) = c("HLA_index",	"HLA_asso_training_size",	"pvalue")
colnames(df_cmv_pval2) = c("HLA_index",	"HLA_asso_training_size",	"pvalue")

df_cmv_pval = rbind(df_cmv_pval1, df_cmv_pval2)
dim(df_cmv_pval)

df_cmv_pval_sorted <- df_cmv_pval[order(df_cmv_pval$HLA_index), ]
stopifnot(all(df_cmv_pval_sorted$HLA_index==(1:220)))

col_cmv_asso_pval = df_cmv_pval_sorted$pvalue

# ---------------------------
# put the table together
# ---------------------------

df_table = data.frame(HLA = col_HLA, 
                      agnostic_AUC = col_HLA_agnostic_AUC,
                      specific_AUC = col_HLA_specific_AUC, 
                      combined_AUC = col_combined_AUC,
                      knn_AUC = col_HLA_weighted_AUC, 
                      specific_training_size = col_HLA_specific_training_size,
                      knn_training_size = col_HLA_weighted_training_size,
                      nTCR_agnostic = col_n_sign_agnostic,
                      nTCR_specific = col_n_sign_specific, 
                      nTCR_knn = col_n_sign_weighted,
                      nTCR_both_agnostic_specific = col_n_intersect_agnostic_specific,
                      nTCR_agnostic_not_specific = col_n_agnostic_only,
                      nTCR_specific_not_agnostic = col_n_specific_only,
                      prop_as_by_a = col_prop_divided_by_agnostic, 
                      prop_as_by_s = col_prop_divided_by_specific,
                      nTCR_both_agnostic_knn = col_n_intersect_agnostic_weighted,
                      nTCR_agnostic_not_knn = col_n_agnostic_not_in_weighted, 
                      nTCR_knn_not_agnostic = col_n_weighted_not_in_agnostic, 
                      prop_ak_by_a = col_both_agnostic_weighted_divided_by_agnostic,
                      prop_ak_by_k = col_both_agnostic_weighted_divided_by_weighted, 
                      HLA_name = col_HLA_name, 
                      HLA_group = col_HLA_group, 
                      HLA_frequency = col_HLA_frequency, 
                      run_10splits = col_run_10splits, 
                      CMV_asso_pvalue=col_cmv_asso_pval)

write.csv(df_table, 
          "../result_tables/full_table.csv", 
          row.names=FALSE)

df_table_kept = df_table[which(df_table$HLA_frequency>=40),]
dim(df_table_kept)

df_table_I = df_table_kept[which(df_table_kept$HLA_group%in%c("A", "B", "C")),]
df_table_II = df_table_kept[which(!(df_table_kept$HLA_group%in%c("A", "B", "C"))),]

df_table_I_sorted = df_table_I[order(df_table_I$HLA), ]
df_table_II_sorted = df_table_II[order(df_table_II$HLA), ]


all_col_names = c("Column Name", 
                  colnames(df_table), 
                  "", 
                  "Note:")

all_col_meanings = c("Interpretation", 
                     "Full name of HLA", 
                     "The AUC on test subjects with the HLA based on prediction of HLA-agnostic model", 
                     "The AUC on test subjects with the HLA based on prediction of HLA-specific model", 
                     "The AUC on test subjects with the HLA based on prediction of combined model",
                     "The AUC on test subjects with the HLA based on prediction of KNN model",
                     "The number of individuals involved in training for HLA-specific model", 
                     "The number of individuals involved in training for KNN model", 
                     "The number of significant TCRs from HLA-agnostic model", 
                     "The number of significant TCRs from HLA-specific model", 
                     "The number of significant TCRs from KNN model", 
                     "The number of TCRs that are significant both from HLA-agnostic model and HLA-specific model", 
                     "The number of TCRs that are significant from HLA-agnostic model, but not from HLA-specific model", 
                     "The number of TCRs that are significant from HLA-specific model, but not from HLA-agnostic model", 
                     "The number of TCRs that are significant both from HLA-agnostic model and HLA-specific model, divided by the number of significant TCRs from HLA-agnostic model", 
                     "The number of TCRs that are significant both from HLA-agnostic model and HLA-specific model, divided by the number of significant TCRs from HLA-specific model", 
                     "The number of TCRs that are significant both from HLA-agnostic model and KNN model", 
                     "The number of TCRs that are significant from HLA-agnostic model, but not from KNN model", 
                     "The number of TCRs that are significant from KNN model, but not from HLA-agnostic model", 
                     "The number of TCRs that are significant both from HLA-agnostic model and KNN model, divided by the number of significant TCRs from HLA-agnostic model", 
                     "The number of TCRs that are significant both from HLA-agnostic model and KNN model, divided by the number of significant TCRs from KNN model", 
                     "The short name of HLA", 
                     "The group of HLA", 
                     "The number of individuals with the HLA, out of all 641 individuals with known CMV status (CMV status not being NA)", 
                     "Whether the HLA goes through 10 training/test splits", 
                     "The association p-value from Fisher's exact test between HLA and CMV")

all_col_meanings = c(all_col_meanings, 
                     "", 
                     paste0("For the HLAs that go through 10 training/test splits (the group of HLAs with frequency between [40, 70]), ", 
                            "their values in terms of models (for example, number of traning individuals, number of significant TCRs, etc.),", 
                            "are the average of metrics based on the splits from which the HLA has at least 32 training and at least 8 test individuals."))

df_readme = data.frame(col1 = all_col_names, 
                       col2 = all_col_meanings)


wb <- createWorkbook()

addWorksheet(wb, "README")
addWorksheet(wb, "HLA I")
addWorksheet(wb, "HLA II")

writeData(wb, sheet = "README", x = df_readme, 
          rowNames=FALSE, colNames = FALSE)
writeData(wb, sheet = "HLA I", x = df_table_I_sorted, 
          rowNames=FALSE, colNames = TRUE)
writeData(wb, sheet = "HLA II", x = df_table_II_sorted, 
          rowNames=FALSE, colNames = TRUE)

# Save the workbook
saveWorkbook(wb, "../result_tables/Table_S1.xlsx", overwrite = TRUE)


sessionInfo()
q(save = "no")

