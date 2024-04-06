# verify whether three versions of HLA similarity matrices are the same

# one is the version under two_gcns folder

df_two_gcns = read.csv(file.path("/Users/sliu/Documents/Fred_Hutch/aTCR_full/two_gcns/data", 
                                 "hla_similarity_blosum62X.csv"), 
                       header=TRUE)
dim(df_two_gcns)

# one is the version under GCN_trial folder

df_gcn_trial = read.csv(file.path("/Users/sliu/Documents/Fred_Hutch/DePTH_folder/GCN_trial/data", 
                                 "hla_similarity_blosum62X.csv"), 
                       header=TRUE)
dim(df_gcn_trial)

summary(c(as.matrix(df_two_gcns - df_gcn_trial)))

# another one is the version under grant folder

list_grants = readRDS(file.path("/Users/sliu/Documents/Fred_Hutch/grants/R01_TCR_revision/results", 
                                "step3_correlation_matrices.rds"))
length(list_grants)
names(list_grants)
list_grants[["aa"]]

summary(c(as.matrix(df_two_gcns) - as.matrix(list_grants[["aa"]])))

# the conclusion is that all three resources are the same
