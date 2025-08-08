
library(readxl)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(ggpointdensity)
library(ggExtra)

theme_set(theme_classic())

fig_dir = "./figures"
res_dir = "../result_tables"

dir_list = list()
dir_list[["hla_i"]] = "../HLA_specific_model_HLA_I"
dir_list[["hla_ii"]] = "../HLA_specific_model_HLA_II"

df_hla_i_index = read.csv("../intermediate_files/HLA_I_index.csv", 
                          header=FALSE)
df_hla_ii_index = read.csv("../intermediate_files/HLA_II_index.csv", 
                           header=FALSE)
df_hla_i_index$hla_type = "hla_i"
df_hla_ii_index$hla_type = "hla_ii"
df_hla_i_index$relative_ind = 0:(nrow(df_hla_i_index)-1)
df_hla_ii_index$relative_ind = 0:(nrow(df_hla_ii_index)-1)

df_hlas = rbind(df_hla_i_index, df_hla_ii_index)
dim(df_hlas)
colnames(df_hlas)[1] = "ind"
head(df_hlas)

df_all_names = read.csv("../intermediate_files/complete_HLA_rownames.csv", 
                        header=TRUE)
dim(df_all_names)
head(df_all_names)

df_high_inds = read.csv("../intermediate_files/high_HLA_index_complete.csv", 
                        header=FALSE)
high_inds = df_high_inds$V1

# record what proportion of CMV associated TCRs under each HLA
# are still significant in terms of TCR-HLA p-value

p_cutoff = 0.001

n_cmv = NULL
n_hla = NULL
n_inter = NULL
# record the spearman correlation between two pvalues
spearman_corrs = NULL 
# record the spearman correlation between two pvalues constrained on CMV pvalue < 0.001
sign_spearman_corrs = NULL 
# record the proportoin of CMV-significant TCRs that are HLA-significant
# under pvalue cutoffs 0.01 and 0.001
prop_mat = NULL 

p_list = list()
cnt = 0

for (i in high_inds){
  
  df_row = df_hlas[which(df_hlas$ind==i),]
  type_i = df_row$hla_type
  relative_i = df_row$relative_ind
  dir_i = dir_list[[type_i]]
  cmv_pval_file = file.path(dir_i, 
                            "results/pvals", 
                            paste0(sprintf("pvalues_hla_index_%d.csv", i)))
  hla_pval_file = file.path(dir_i, 
                            "results/TCR_HLA_pvalues", 
                            paste0(as.character(relative_i), "_TCR_HLA_pvalues.csv")) 
  df_cmv_pval = read.csv(cmv_pval_file, header=TRUE)
  df_hla_pval = read.csv(hla_pval_file, header=TRUE)  
  
  cmv_sign_tcrs = which(df_cmv_pval$pval<p_cutoff) 
  
  n_cmv = c(n_cmv, length(cmv_sign_tcrs))
  n_hla = c(n_hla, sum(df_hla_pval$pval < p_cutoff))
  n_inter = c(n_inter, 
              sum((df_cmv_pval$pval < p_cutoff) & (df_hla_pval$pval < p_cutoff)))
  
  df_plot = data.frame(cmv_pval = df_cmv_pval$pval, 
                       hla_pval = df_hla_pval$pval)
  
  spearman_corrs = c(spearman_corrs, 
                     cor(df_cmv_pval$pval, 
                         df_hla_pval$pval, 
                         method = "spearman"))
    
  df_plot = df_plot[cmv_sign_tcrs, ]

  sign_spearman_corrs = c(sign_spearman_corrs, 
                          cor(df_plot$cmv_pval, 
                              df_plot$hla_pval, 
                              method = "spearman"))
  
  cur_vec = c(mean(df_plot$hla_pval<0.01), 
              mean(df_plot$hla_pval<0.001))
  
  prop_mat = rbind(prop_mat, cur_vec)
  
  cnt = cnt + 1
  
  p_list[[cnt]] = ggplot(df_plot, aes(x = -log10(cmv_pval), y = -log10(hla_pval))) +
                  geom_pointdensity() +
                  scale_color_viridis() + 
                  geom_hline(yintercept=-log10(p_cutoff), color = "grey", linetype = "dashed") + 
                  ggtitle(paste0(df_all_names$HLA_name[i+1], "\neach point is one TCR with CMV pval <0.001", 
                                 "\n", as.character(round(100*n_inter[length(n_inter)]/n_cmv[length(n_cmv)], digits=1)), 
                                 "% TCRs sign with CMV still sign with HLA"))
  
}

figure_file = file.path(fig_dir, 
                        paste0("TCR_HLA_CMV_pvalue_scatterplots_high_indexes.pdf"))

n_col = 4
n_row = ceiling(length(p_list)/n_col)

pdf(file = figure_file, 
    width = 4.5*n_col, height = 3.5*n_row+1)
combined_plot = ggarrange(plotlist = p_list, ncol = n_col, nrow = n_row)
final_plot <- annotate_figure(combined_plot,
                              top = text_grob(paste0("Scatterplots between -log10(pvals) with CMV and with HLA\nconstrained on TCR having pval with CMV < 0.001"), 
                                              size = 14, face = "bold"))
print(final_plot)
dev.off()    

n_inter/n_cmv

df_high = data.frame(n_sign_w_CMV = n_cmv, 
                     n_sign_w_HLA = n_hla, 
                     n_sign_w_both = n_inter,
                     prop_out_of_CMV = n_inter/n_cmv, 
                     prop_out_of_HLA = n_inter/n_hla)

# load training individual indexes
df_train_index = read.csv("../intermediate_files/train_index.csv", 
                          header=FALSE)
min(df_train_index$V1)

# load full HLA matrix with split HLA-II haplotypes and NAs kepts
HLA_matrix = read.csv("../intermediate_files/complete_HLA_NAs_kept.csv")
HLA_train_subset = HLA_matrix[, (df_train_index$V1+1)]

HLA_known = rowSums(!is.na(HLA_train_subset))
HLA_pos = rowSums(HLA_train_subset==1, na.rm=TRUE)

table(HLA_known[high_inds+1])
summary(HLA_pos[high_inds+1])

df_high$n_persons_known = HLA_known[high_inds+1]
df_high$n_persons_pos = HLA_pos[high_inds+1]

q_list = list()

q_list[[1]] = ggplot(df_high, aes(x = log10(n_persons_known), y = prop_out_of_CMV)) +
  geom_pointdensity() +
  scale_color_viridis() + 
  ggtitle(paste0("prop of TCR sign with CMV still sign with HLA", 
                 "\nagainst log10(train subjects \nwith status being not NA for HLA)"))

q_list[[2]] = ggplot(df_high, aes(x = log10(n_persons_pos), y = prop_out_of_CMV)) +
  geom_pointdensity() +
  scale_color_viridis() + 
  ggtitle(paste0("prop of TCR sign with CMV still sign with HLA", 
                 "\nagainst log10(train subjects \nhaving the HLA)"))

q_list[[3]] = ggplot(df_high, aes(x = log10(n_persons_known), y = log10(n_persons_pos))) +
  geom_pointdensity() +
  scale_color_viridis() + 
  ggtitle(paste0("log10(train subjects having the HLA)", 
                 "\nagainst log10(train subjects \nwith status being not NA for HLA)"))

figure_file = file.path(fig_dir, 
                        paste0("prop_of_TCR_sign_with_CMV_kept_high_indexes.pdf"))

n_col = 2
n_row = ceiling(length(q_list)/n_col)

pdf(file = figure_file, 
    width = 4.5*n_col, height = 3.5*n_row)
combined_plot = ggarrange(plotlist = q_list, ncol = n_col, nrow = n_row)
print(combined_plot)
dev.off()    

# out put these items:
# (1) the correlation between the CMV pvalues and HLA pvalues
# (2) the proportion of the CMV-significant TCRs (pvalue < 0.001) 
#     that are associated with the HLA at different pvalue cutoffs (0.01, 0.001)

summary(spearman_corrs)
summary(sign_spearman_corrs)

row.names(prop_mat) = df_all_names$HLA_name[high_inds+1]
colnames(prop_mat) = c("prop_pvals_below_0.01", "prop_pvals_below_0.001")
write.csv(prop_mat, 
          file = file.path(res_dir, 
                           "Table_S3.csv"), 
          row.names=TRUE)

sessionInfo()
q(save = "no")
