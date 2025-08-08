
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggpointdensity)
library(ggExtra)

options(bitmapType = "cairo") 

theme_set(theme_classic())

# generate scatterplots between the pvalues from A0201 only individuals and 
# those from individuals with both A0201 and B0702

df_hla_names = read.csv("../intermediate_files/complete_HLA_rownames.csv", 
                        header=TRUE)
dim(df_hla_names)

two_hlas = c("HLA-A*02:01", "HLA-B*07:02")
two_hla_indexes = match(two_hlas, df_hla_names$HLA_name)
table(two_hlas==df_hla_names$HLA_name[two_hla_indexes], 
      useNA="ifany")

hla_index_dict = list()
for (i in 1:2){
  hla_index_dict[[two_hlas[i]]] = two_hla_indexes[i]
}

hla_index_dict

# load pvalue files


hla_name = two_hlas[1]

df_pval_A0201 = read.csv(file.path("../HLA_specific_model_HLA_I/results/pvals", 
                                  paste0("pvalues_hla_index_", 
                                         as.character(hla_index_dict[[hla_name]]-1), ".csv")), 
                         header = TRUE)
stopifnot(dim(df_pval_A0201)[1]==1098738)
head(df_pval_A0201)

df_pval_double = read.csv("./results/pvalues_double_hlas.csv", 
                          header = TRUE)
stopifnot(dim(df_pval_double)[1]==1098738)
head(df_pval_double)



df_cur = data.frame(A0201=df_pval_A0201$pval, 
                    double=df_pval_double$pval)

# df_cur = df_cur[which((df_cur$A0201<1)|(df_cur$double<1)),]
  
p1 <- ggplot(df_cur, aes(x = -log10(A0201), y = -log10(double))) +
      geom_pointdensity() +
      scale_color_viridis() + 
      xlab("-log10(A0201-specific pvalue)") + 
      ylab(paste0("-log10(A0201-B0702-specific pvalue)")) +     
      geom_abline(intercept = 0, slope = 1,  color = "grey", linetype = "dashed") +
      ggtitle(paste0("comparison between pvalues from \n", 
                     "A0201-specific and A0201-B0702-specific model\n", 
                     "each point is one TCR"))
  

figure_dir = "./figures"
if (!dir.exists(figure_dir)){
  dir.create(figure_dir)
}

figure_file = file.path(figure_dir, 
                        paste0("TCR_pval_single_double_models_scatterplot.png"))

png(file = figure_file, 
    width = 6.8, height = 5.8, units="in", res=1200)
print(p1)
dev.off()    


sessionInfo()
q(save = "no")