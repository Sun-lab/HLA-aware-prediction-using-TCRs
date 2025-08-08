# compare the new and old table S1 in terms of overlapped HLAs

library("readxl")
library("ggplot2")

new_HLA_I <- read_excel("../result_tables/Table_S1.xlsx", sheet = "HLA I")
new_HLA_II <- read_excel("../result_tables/Table_S1.xlsx", sheet = "HLA II")

old_HLA_I <- read_excel("../old_table/Table_S1.xlsx", sheet = "HLA I")
old_HLA_II <- read_excel("../old_table/Table_S1.xlsx", sheet = "HLA II")

dim(old_HLA_I)
dim(old_HLA_II)

old_HLA_I_matched = old_HLA_I[match(new_HLA_I$HLA, old_HLA_I$HLA),]
old_HLA_II_matched = old_HLA_II[match(new_HLA_II$HLA, old_HLA_II$HLA),]

summary(new_HLA_I$specific_AUC-old_HLA_I_matched$HLA_specific_AUC)
summary(new_HLA_I$agnostic_AUC-old_HLA_I_matched$HLA_ignorant_AUC)
summary(new_HLA_I$knn_AUC-old_HLA_I_matched$KNN_AUC)

summary(new_HLA_II$specific_AUC-old_HLA_II_matched$HLA_specific_AUC)
table(new_HLA_II$specific_AUC!=old_HLA_II_matched$HLA_specific_AUC)

summary(new_HLA_II$agnostic_AUC-old_HLA_II_matched$HLA_ignorant_AUC)
table(new_HLA_II$agnostic_AUC!=old_HLA_II_matched$HLA_ignorant_AUC)

summary(new_HLA_II$knn_AUC-old_HLA_II_matched$KNN_AUC)
table(new_HLA_II$knn_AUC!=old_HLA_II_matched$KNN_AUC)

new_HLA_I_false = new_HLA_I[new_HLA_I$run_10splits==FALSE,]
new_HLA_II_false = new_HLA_II[new_HLA_II$run_10splits==FALSE,]

old_HLA_I_matched_false = old_HLA_I_matched[new_HLA_I$run_10splits==FALSE,]
old_HLA_II_matched_false = old_HLA_II_matched[new_HLA_II$run_10splits==FALSE,]

summary(abs(new_HLA_I_false$agnostic_AUC - old_HLA_I_matched_false$HLA_ignorant_AUC))
summary(abs(new_HLA_I_false$specific_AUC - old_HLA_I_matched_false$HLA_specific_AUC))
summary(abs(new_HLA_I_false$specific_training_size - old_HLA_I_matched_false$Training_Size))
summary(abs(new_HLA_I_false$knn_AUC - old_HLA_I_matched_false$KNN_AUC))
summary(abs(new_HLA_I_false$knn_training_size - old_HLA_I_matched_false$KNN_training_size))

summary(abs(new_HLA_II_false$agnostic_AUC - old_HLA_II_matched_false$HLA_ignorant_AUC))
summary(abs(new_HLA_II_false$specific_AUC - old_HLA_II_matched_false$HLA_specific_AUC))
summary(abs(new_HLA_II_false$specific_training_size - old_HLA_II_matched_false$Training_Size))
summary(abs(new_HLA_II_false$knn_AUC - old_HLA_II_matched_false$KNN_AUC))
summary(abs(new_HLA_II_false$knn_training_size - old_HLA_II_matched_false$KNN_training_size))

ready_HLA_I <- read.csv("../../plot/HLA_I_ready_data.csv", header=TRUE)
ready_HLA_II <- read.csv("../../plot/HLA_II_ready_data.csv", header=TRUE)

dim(ready_HLA_I)
dim(ready_HLA_II)

ready_HLA_I_matched = ready_HLA_I[match(new_HLA_I$HLA, ready_HLA_I$HLA),]
ready_HLA_II_matched = ready_HLA_II[match(new_HLA_II$HLA, ready_HLA_II$HLA),]

ready_HLA_I_matched = ready_HLA_I[match(new_HLA_I$HLA, ready_HLA_I$HLA),]
ready_HLA_II_matched = ready_HLA_II[match(new_HLA_II$HLA, ready_HLA_II$HLA),]

ready_HLA_I_matched_false = ready_HLA_I_matched[new_HLA_I$run_10splits==FALSE,]
ready_HLA_II_matched_false = ready_HLA_II_matched[new_HLA_II$run_10splits==FALSE,]

table(new_HLA_I_false$HLA==ready_HLA_I_matched_false$HLA, useNA="ifany")
table(new_HLA_II_false$HLA==ready_HLA_II_matched_false$HLA, useNA="ifany")

summary(abs(new_HLA_I_false$combined_AUC - ready_HLA_I_matched_false$Combined_AUC))
summary(abs(new_HLA_II_false$combined_AUC - ready_HLA_II_matched_false$Combined_AUC))

ggplot(new_HLA_I_false, aes(x=specific_training_size, y=combined_AUC-agnostic_AUC)) + 
  geom_point() + 
  xlab("training sample size") + 
  ylab("combined - agnostic") +
  theme_classic() + ggtitle("HLA-I")

ggplot(new_HLA_II_false, aes(x=specific_training_size, y=combined_AUC-agnostic_AUC)) + 
  geom_point() + 
  xlab("training sample size") + 
  ylab("combined - agnostic") +
  theme_classic() + ggtitle("HLA-II")


df_ready_I = read.csv("../../plot/HLA_I_ready_data.csv", header=TRUE)
df_ready_II = read.csv("../../plot/HLA_II_ready_data.csv", header=TRUE)

df_ready_I_matched = df_ready_I[match(new_HLA_I$HLA, df_ready_I$HLA),]
df_ready_II_matched = df_ready_II[match(new_HLA_II$HLA, df_ready_II$HLA),]

table(new_HLA_I$HLA==df_ready_I_matched$HLA, useNA="ifany")
table(new_HLA_II$HLA==df_ready_II_matched$HLA, useNA="ifany")

new_HLA_I$old_combined = df_ready_I_matched$Combined_AUC
new_HLA_II$old_combined = df_ready_II_matched$Combined_AUC

new_HLA_I_false = new_HLA_I[new_HLA_I$run_10splits==FALSE,]
new_HLA_II_false = new_HLA_II[new_HLA_II$run_10splits==FALSE,]

summary(new_HLA_I_false$combined_AUC - new_HLA_I_false$old_combined)
summary(new_HLA_II_false$combined_AUC - new_HLA_II_false$old_combined)

sessionInfo()
q(save = "no")

