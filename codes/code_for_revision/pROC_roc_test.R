# run pROC/roc.test for HLAs with population frequency greater than 70

library(pROC)

df_high = read.csv("../intermediate_files/high_HLA_index_complete.csv", 
                   header=FALSE)
dim(df_high)
head(df_high)

high_ids = df_high$V1

df_hla_names = read.csv("../intermediate_files/complete_HLA_rownames.csv", 
                        header=TRUE)
hlas = df_hla_names$HLA_name

auc_combined = rep(NA, nrow(df_high))
auc_agnostic = rep(NA, nrow(df_high))
DeLong_pvalues = rep(NA, nrow(df_high))

cnt = 0

for (i in high_ids){
  
  print(sprintf("HLA index %d", i))
  
  cnt = cnt + 1
  
  df_roc_combined = read.csv(file.path("../combined_HLA_model/results/roc_related", 
                                       sprintf("roc_related_%d.csv", i)), header=TRUE)
  roc_combined = roc(df_roc_combined$label, df_roc_combined$prediction)
  auc_combined[cnt] = auc(roc_combined)
  
  df_roc_agnostic = read.csv(file.path("../HLA_agnostic_model/results/roc_related", 
                                       sprintf("roc_related_%d.csv", i)), header=TRUE)
  roc_agnostic = roc(df_roc_agnostic$label, df_roc_agnostic$prediction)    
  auc_agnostic[cnt] = auc(roc_agnostic)
  
  roc.test(roc_combined, roc_agnostic, alternative="greater")$p.value
  

  DeLong_pvalues[cnt] = roc.test(roc_combined, roc_agnostic, alternative="greater")$p.value
  
}

df = data.frame(hla_name = hlas[high_ids+1], 
                auc_combined = auc_combined, 
                auc_agnostic = auc_agnostic, 
                DeLong_pvalue = DeLong_pvalues)
                
df_table = read.csv("../result_tables/full_table.csv", header=TRUE)

df_table_matched = df_table[match(df$hla_name, df_table$HLA), ]

stopifnot(all(df$hla_name==df_table_matched$HLA))

min(df_table_matched$HLA_frequency)

summary(df$auc_combined - df_table_matched$combined_AUC)
summary(df$auc_agnostic - df_table_matched$agnostic_AUC)

write.csv(df, 
          file = "results/pROC_roc_test_results.csv", 
          row.names=FALSE)


sessionInfo()
q(save = "no")



