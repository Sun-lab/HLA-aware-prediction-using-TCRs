library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(gridExtra)
library(ggrepel)
library(grid)


get_HLA_II_asso = function()
{
  HLA_II_index = read.csv("../intermediate_files/HLA_II_index.csv",header = F)$V1 + 1
  HLA_II_asso_res = numeric(length = length(HLA_II_index))
  HLA_II_asso_training_size = numeric(length=length(HLA_II_index))
  
  CMV_subj = which(!is.na( CMV ))
  res_index = 1
  for (HLA_index in HLA_II_index) {
    # the subject index for fisher exact test
    subj_inx = intersect(which(!is.na(split_HLA[HLA_index,])), CMV_subj)
    subj_size = length(subj_inx)
    HLA_dat = split_HLA[HLA_index, subj_inx]
    CMV_dat = CMV[subj_inx]
    res_00 = res_01 = res_10 = res_11 = 0;
    
    for (i in 1:subj_size) {
      if (HLA_dat[i] == 0 &&  CMV_dat[i] == 0 ) {
        res_00 = res_00 + 1
      }
      if (HLA_dat[i] == 0 &&  CMV_dat[i] == 1 ) {
        res_01 = res_01 + 1
      }
      if (HLA_dat[i] == 1 &&  CMV_dat[i] == 0 ) {
        res_10 = res_10 + 1
      }
      if (HLA_dat[i] == 1 &&  CMV_dat[i] == 1 ) {
        res_11 = res_11 + 1
      }
    }
    
    test_mat = matrix(c(res_00, res_01, res_10,res_11),nrow=2, byrow = T)
    test_res = fisher.test(test_mat, alternative = "two.sided")
    HLA_II_asso_res[res_index] = test_res$p.value
    HLA_II_asso_training_size[res_index] = subj_size
    res_index = res_index + 1
  }
  # here the indexes are already converted to 1-indexed  
  df_HLA_II_res = data.frame(HLA_II_index = HLA_II_index, 
                             HLA_II_asso_training_size = HLA_II_asso_training_size, 
                             pvalue = HLA_II_asso_res)  
  return(df_HLA_II_res)
}

get_HLA_I_asso = function()
{
  
  HLA_I_index = read.csv("../intermediate_files/HLA_I_index.csv",header = F)$V1 + 1
  HLA_I_asso_res = numeric(85)
  HLA_I_asso_training_size = numeric(85)
  
  CMV_subj = which(!is.na( CMV ))
  res_index = 1
  
  for (HLA_index in HLA_I_index) {
    # the subject index for fisher exact test
    subj_inx = intersect(which(!is.na(split_HLA[HLA_index,])), CMV_subj)
    subj_size = length(subj_inx)
    HLA_dat = split_HLA[HLA_index, subj_inx]
    CMV_dat = CMV[subj_inx]
    res_00 = res_01 = res_10 = res_11 = 0;
    
    for (i in 1:subj_size) {
      if (HLA_dat[i] == 0 &&  CMV_dat[i] == 0 ) {
        res_00 = res_00 + 1
      }
      if (HLA_dat[i] == 0 &&  CMV_dat[i] == 1 ) {
        res_01 = res_01 + 1
      }
      if (HLA_dat[i] == 1 &&  CMV_dat[i] == 0 ) {
        res_10 = res_10 + 1
      }
      if (HLA_dat[i] == 1 &&  CMV_dat[i] == 1 ) {
        res_11 = res_11 + 1
      }
    }
    
    test_mat = matrix(c(res_00, res_01, res_10,res_11),nrow=2, byrow = T)
    test_res = fisher.test(test_mat, alternative = "two.sided")
    HLA_I_asso_res[res_index] = test_res$p.value
    HLA_I_asso_training_size[res_index] = subj_size
    res_index = res_index + 1
  }
  # here the indexes are already converted to 1-indexed
  df_HLA_I_res = data.frame(HLA_I_index = HLA_I_index, 
                            HLA_I_asso_training_size = HLA_I_asso_training_size, 
                            pvalue = HLA_I_asso_res)
  
  return(df_HLA_I_res)
}


split_HLA = readRDS("../../../intermediate_files/complete_HLA_NAs_kept.rds")
CMV = read.csv("../../../specific_HLA_pred/data_files/CMV.txt")
CMV = CMV$x

df_HLA_I_res = get_HLA_I_asso()
df_HLA_II_res = get_HLA_II_asso()


write.csv(df_HLA_I_res, 
          "results/HLA_I_CMV_asso_res.csv", 
          row.names = FALSE)

write.csv(df_HLA_II_res, 
          "results/HLA_II_CMV_asso_res.csv", 
          row.names = FALSE)

sessionInfo()
q(save = "no")
