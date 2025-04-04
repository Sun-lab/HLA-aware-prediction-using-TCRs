
# compute the sequence score matrix based on blosum62 using python
# load the results and compute distance based on the formula in 
# NetMHCpan paper

HLA_filepath = "../../data/HLA_I_pseudo_40.csv"
df_HLA = read.csv(HLA_filepath, header = TRUE)

df_blosum62_X_score = 
  read.csv("../../data/blosum62_X_seq_align_score_HLA_I.csv", header = TRUE)
df_blosum62_X_score[1:6, 1:6]

# the column order of HLA matches that from HLA pseudo sequence file
simi_mat_blosum62_X = matrix(NA, ncol = nrow(df_HLA)[1], nrow = nrow(df_HLA)[1])

for (i in 1:nrow(df_blosum62_X_score)){
  for (j in 1:nrow(df_blosum62_X_score)){
    simi_mat_blosum62_X[i, j] = 
      df_blosum62_X_score[i,j]/sqrt(df_blosum62_X_score[i,i]*df_blosum62_X_score[j,j])
  }
}

rownames(simi_mat_blosum62_X) = df_HLA$hla

corr_matrices = list()
corr_matrices[["aa"]] = simi_mat_blosum62_X

saveRDS(corr_matrices, 
        file = "../../data/correlation_matrices.rds")


sessionInfo()
q(save="no")
