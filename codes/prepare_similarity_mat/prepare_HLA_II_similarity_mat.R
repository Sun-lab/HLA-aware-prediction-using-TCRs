

HLA_filepath = "../../data/HLA_II_pseudo_45.csv"
df_HLA = read.csv(HLA_filepath, header = TRUE)
dim(df_HLA)
head(df_HLA)

df_blosum62_X_score = 
  read.csv("../../data/blosum62_X_seq_align_score_HLA_II.csv", header = TRUE)
dim(df_blosum62_X_score)    
df_blosum62_X_score[1:6, 1:6]

# the column order of HLA matches that from HLA pseudo sequence file
simi_mat_blosum62_X = matrix(NA, ncol = nrow(df_HLA)[1], nrow = nrow(df_HLA)[1])

for (i in 1:nrow(df_blosum62_X_score)){
  for (j in 1:nrow(df_blosum62_X_score)){
    simi_mat_blosum62_X[i, j] = 
      df_blosum62_X_score[i,j]/sqrt(df_blosum62_X_score[i,i]*df_blosum62_X_score[j,j])
  }
}

summary(c(simi_mat_blosum62_X))

rownames(simi_mat_blosum62_X) = df_HLA$hla

saveRDS(simi_mat_blosum62_X, 
        file = "../../data/HLA_II_similarity_matrix_aa.rds")


sessionInfo()
q(save="no")
