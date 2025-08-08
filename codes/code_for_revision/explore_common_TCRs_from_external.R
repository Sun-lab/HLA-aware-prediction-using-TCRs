
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggpointdensity)
library(ggExtra)

theme_set(theme_classic())

# for external TCRs from Chen et al. 2017 and Huth et al. 2019 that also appear
# in the Emerson frequent TCRs, 
# identify their index among the Emerson frequent TCRs
# load the pvalue files from HLA-agnostic model and HLA-specific models under corresponding HLAs
# visualize the pvalue scatterplots


# load external TCRs from Chen et al. 2017 and Huth et al. 2019 that also appear
# in the Emerson frequent TCRs

df_A02 = read.csv(file.path("../intermediate_files/Chen_2017", 
                            "Chen_2017_TCRs_processed_common.csv"), 
                  header=TRUE)
dim(df_A02)
head(df_A02)

df_B0702 = read.csv(file.path("../intermediate_files/Huth_2019", 
                            "Huth_2019_TCRs_processed_common_B0702.csv"), 
                  header=TRUE)
dim(df_B0702)

df_C0702 = read.csv(file.path("../intermediate_files/Huth_2019", 
                              "Huth_2019_TCRs_processed_common_C0702.csv"), 
                    header=TRUE)
dim(df_C0702)

# load Emerson TCR names

df_names = read.csv(file.path("../intermediate_files", 
                              "TCR_names.csv"), 
                    header=TRUE)
dim(df_names)
head(df_names)
  
# find index
A02_indexes = match(df_A02$TCR_name, df_names$TCR_name)
B0702_indexes = match(df_B0702$TCR_name, df_names$TCR_name)
C0702_indexes = match(df_C0702$TCR_name, df_names$TCR_name)

table(df_names$TCR_name[A02_indexes]==df_A02$TCR_name, useNA="ifany")
table(df_names$TCR_name[B0702_indexes]==df_B0702$TCR_name, useNA="ifany")
table(df_names$TCR_name[C0702_indexes]==df_C0702$TCR_name, useNA="ifany")

# save 0-indexed indexes out, for the consideration of running python
output_A02 = data.frame(TCR_index = A02_indexes-1)
write.csv(output_A02, 
          file.path("../intermediate_files/Chen_2017", 
                    "Chen_2017_TCRs_common_indexes_A02.csv"),
          row.names=FALSE)

output_B0702 = data.frame(TCR_index = B0702_indexes-1)
write.csv(output_B0702, 
          file.path("../intermediate_files/Huth_2019", 
                    "Huth_2019_TCRs_common_indexes_B0702.csv"),
          row.names=FALSE)

output_C0702 = data.frame(TCR_index = C0702_indexes-1)
write.csv(output_C0702, 
          file.path("../intermediate_files/Huth_2019", 
                    "Huth_2019_TCRs_common_indexes_C0702.csv"),
          row.names=FALSE)

# identify the index of three HLAs, to help with loading pvalue files

df_hla_names = read.csv("../intermediate_files/complete_HLA_rownames.csv", 
                        header=TRUE)
dim(df_hla_names)

three_hlas = c("HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02")
three_hla_indexes = match(three_hlas, df_hla_names$HLA_name)
table(three_hlas==df_hla_names$HLA_name[three_hla_indexes], 
      useNA="ifany")

hla_index_dict = list()
for (i in 1:3){
  hla_index_dict[[three_hlas[i]]] = three_hla_indexes[i]
}

# load pvalue files

pvalue_list = list()

df_pval_all = read.csv("../HLA_agnostic_model/results/pvalues.csv", 
                       header = TRUE)
stopifnot(dim(df_pval_all)[1]==1098738)
pvalue_list[["all"]] = df_pval_all$pval

for (hla_name in three_hlas){
  df_pval_specific = read.csv(file.path("../HLA_specific_model_HLA_I/results/pvals", 
                                        paste0("pvalues_hla_index_", 
                                               as.character(hla_index_dict[[hla_name]]-1), ".csv")), 
                              header = TRUE)
  stopifnot(dim(df_pval_specific)[1]==1098738)
  pvalue_list[[hla_name]] = df_pval_specific$pval
}


TCR_index_list = list()

TCR_index_list[["HLA-A*02:01"]] = A02_indexes
TCR_index_list[["HLA-B*07:02"]] = B0702_indexes
TCR_index_list[["HLA-C*07:02"]] = C0702_indexes

fig_list = list()

cnt = 0

for (hla_name in three_hlas){

  cur_indexes = TCR_index_list[[hla_name]]
  df_cur = data.frame(all=pvalue_list[["all"]][cur_indexes], 
                      hla=pvalue_list[[hla_name]][cur_indexes])
  
    
  p1 <- ggplot(df_cur, aes(x = -log10(all), y = -log10(hla))) +
        geom_pointdensity() +
        scale_color_viridis() + 
        xlab("-log10(HLA-agnostic pvalue)") + 
        ylab(paste0("-log10(", hla_name, "-specific pvalue)")) +     
        geom_abline(intercept = 0, slope = 1,  color = "grey", linetype = "dashed") +
        ggtitle(paste0(hla_name, " TCRs\nappearing in Emerson"))
  
  cnt = cnt + 1
  fig_list[[cnt]] = p1
  
}

# also add the scatterplots based on the TCRs significant under each HLA-specific model

sign_TCR_list = list()

for (hla_name in three_hlas){
  sign_TCR_list[[hla_name]] = which(pvalue_list[[hla_name]]<=0.001)
}


for (hla_name in three_hlas){
  
  cur_indexes = sign_TCR_list[[hla_name]]
  df_cur = data.frame(all=pvalue_list[["all"]][cur_indexes], 
                      hla=pvalue_list[[hla_name]][cur_indexes])
  
  
  p1 <- ggplot(df_cur, aes(x = -log10(all), y = -log10(hla))) +
    geom_pointdensity() +
    scale_color_viridis() + 
    xlab("-log10(HLA-agnostic pvalue)") + 
    ylab(paste0("-log10(", hla_name, "-specific pvalue)")) +     
    geom_abline(intercept = 0, slope = 1,  color = "grey", linetype = "dashed") +
    ggtitle(paste0(hla_name, "\nHLA-specific model significant TCRs"))
  
  cnt = cnt + 1
  fig_list[[cnt]] = p1
  
}

# get from the TCRs significant under each HLA-specific model
# under each pvalue cutoffs
# what proportion of them are also in the external TCRs

sign_TCR_cutoffs = list()

p_cutoffs = 10^c(-3, -4, -5, -6)

for (pcut in p_cutoffs){
  sign_TCR_list = list()
  for (hla_name in three_hlas){
    sign_TCR_list[[hla_name]] = which(pvalue_list[[hla_name]]<=pcut)
  }
  sign_TCR_cutoffs[[as.character(pcut)]] = sign_TCR_list
}

for (hla_name in three_hlas){
  
  prop_vec = NULL
  
  for (pcut in p_cutoffs){
    
    cur_sign_TCRs = sign_TCR_cutoffs[[as.character(pcut)]][[hla_name]]
    cur_intersect = intersect(cur_sign_TCRs, TCR_index_list[[hla_name]])
  
    prop_vec = c(prop_vec, length(cur_intersect)/length(cur_sign_TCRs))
    
  }
  
  df_prop = data.frame(p_cutoffs = c("1e-03", "1e-04", "1e-05", "1e-06"), 
                       prop = prop_vec)
    
  p1 <- ggplot(df_prop, aes(fill=p_cutoffs, y=prop, x=p_cutoffs)) + 
        geom_bar(position="dodge", stat="identity") + 
        xlab("p-value cutoff") + 
        ylab("Proportion") + 
        ylim(0, 1) + 
        ggtitle(paste0(hla_name, "\nProportion of external TCRs out of", 
                       "\nTCRs significant from HLA-specific model"))
  
  cnt = cnt + 1
  fig_list[[cnt]] = p1
  
}

figure_dir = "./figures"
dir.create(figure_dir)

figure_file = file.path(figure_dir, 
                        paste0("external_TCRs_pvalue_scatterplots.pdf"))

n_col = 3
n_row = ceiling(length(fig_list)/n_col)

pdf(file = figure_file, 
    width = 4.2*n_col, height = 3.2*n_row)
combined_plot = ggarrange(plotlist = fig_list, ncol = n_col, nrow = n_row)
print(combined_plot)
dev.off()    

# get what proportion of external TCRs appearing in Emerson data appear significant
# under HLA-specific model

prop_sign_vec = NULL

for (hla_name in three_hlas){
  
  pcut = p_cutoffs[1]
  
  cur_sign_TCRs = sign_TCR_cutoffs[[as.character(pcut)]][[hla_name]]
  cur_intersect = intersect(cur_sign_TCRs, TCR_index_list[[hla_name]])
  
  prop_sign_vec = c(prop_sign_vec, length(cur_intersect)/length(TCR_index_list[[hla_name]]))
    
}

prop_sign_vec 


sessionInfo()
q(save = "no")