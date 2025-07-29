
library(readxl)

df_table_i = read_excel("../result_tables/Table_S1.xlsx", sheet = "HLA I")
df_table_i = as.data.frame(df_table_i)
head(df_table_i)

df_table_ii = read_excel("../result_tables/Table_S1.xlsx", sheet = "HLA II")
df_table_ii = as.data.frame(df_table_ii)
head(df_table_ii)

df_hla_i_names = read.csv("../intermediate_files/HLA_I_names.csv", header=TRUE)
df_hla_ii_names = read.csv("../intermediate_files/HLA_II_names.csv", header=TRUE)

df_hla_i_index = read.csv("../intermediate_files/HLA_I_index.csv", header=FALSE)
df_hla_ii_index = read.csv("../intermediate_files/HLA_II_index.csv", header=FALSE)

hla_i_dir = "../HLA_specific_model_HLA_I/results/permutation_aucs"
hla_ii_dir = "../HLA_specific_model_HLA_II/results/permutation_aucs"

six_hla_names = NULL
six_pvalues = NULL

hla_i_inds = c(23)

for (hla_i in hla_i_inds){
  hla_name = df_hla_i_names$x[hla_i+1]
  real_auc = df_table_i$specific_AUC[which(df_table_i$HLA==hla_name)]
  
  permute_aucs = NULL
  for (j in 0:99){
    df_j  = read.csv(file.path(hla_i_dir, 
                               paste0("hla_index_", as.character(hla_i)), 
                               paste0(as.character(hla_i), "_", as.character(j), "_chunk_res.csv")), 
                     header=TRUE)
    permute_aucs = c(permute_aucs, df_j$permutation_auc)
  }
  print(real_auc)
  print(summary(permute_aucs))
  stopifnot(sum(permute_aucs==-1)==0)
  pvalue_auc = mean(permute_aucs>=real_auc)
  print(sum(is.na(permute_aucs)))
  print(paste0(hla_name, " permutation p-value: ", as.character(pvalue_auc)))
  
  six_hla_names = c(six_hla_names, hla_name)
  six_pvalues = c(six_pvalues, pvalue_auc)
  
}


hla_ii_inds = c(21, 28, 33, 65, 130)

for (hla_ii in hla_ii_inds){
  hla_name = df_hla_ii_names$x[hla_ii+1]
  real_auc = df_table_ii$specific_AUC[which(df_table_ii$HLA==hla_name)]
  
  permute_aucs = NULL
  for (j in 0:99){
    df_j  = read.csv(file.path(hla_ii_dir, 
                               paste0("hla_index_", as.character(hla_ii)), 
                               paste0(as.character(hla_ii), "_", as.character(j), "_chunk_res.csv")), 
                     header=TRUE)
    permute_aucs = c(permute_aucs, df_j$permutation_auc)
  }
  print(real_auc)
  print(summary(permute_aucs))
  stopifnot(sum(permute_aucs==-1)==0)
  pvalue_auc = mean(permute_aucs>=real_auc)
  print(sum(is.na(permute_aucs)))
  print(paste0(hla_name, " permutation p-value: ", as.character(pvalue_auc)))

  six_hla_names = c(six_hla_names, hla_name)
  six_pvalues = c(six_pvalues, pvalue_auc)
  
}

df_output = data.frame(hla_name = six_hla_names, 
                       permutation_pvalue = six_pvalues)

write.csv(df_output, 
          file = "../result_tables/Table_S2.csv", 
          row.names=FALSE)


sessionInfo()
q(save = "no")
