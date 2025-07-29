library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)
library(quantreg)
library(readxl)


get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

df1 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA I")
df2 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA II")

baseline = 0.9163461538461538
colors_ <- c("#1f77b4", "#ff7f0e", "#2ca02c")

# load the HLA-CMV pvalues
df_hla_i_res = read.csv("../../codes/HLA_CMV_association/results/HLA_I_CMV_asso_res.csv", 
                        header=TRUE)
df_hla_ii_res = read.csv("../../codes/HLA_CMV_association/results/HLA_II_CMV_asso_res.csv", 
                         header=TRUE)

df_hla_names = read.csv("../../codes/intermediate_files/complete_HLA_rownames.csv",
                        header=TRUE)
hla_names = df_hla_names$HLA_name

df_hla_i_res$HLA_name = hla_names[df_hla_i_res$HLA_I_index] 
df_hla_ii_res$HLA_name = hla_names[df_hla_ii_res$HLA_II_index] 

df_hla_i_res_matched = df_hla_i_res[match(df1$HLA, df_hla_i_res$HLA_name), ]
df_hla_ii_res_matched = df_hla_ii_res[match(df2$HLA, df_hla_ii_res$HLA_name), ]

table(df1$HLA==df_hla_i_res_matched$HLA_name, useNA="ifany")
table(df2$HLA==df_hla_ii_res_matched$HLA_name, useNA="ifany")

df1$HLA_CMV_asso_pval = df_hla_i_res_matched$pvalue
df2$HLA_CMV_asso_pval = df_hla_ii_res_matched$pvalue

########################
# Specific-AUC v.s. pval
########################

p1 = ggplot(df1, aes(x = -log10(HLA_CMV_asso_pval) , y = specific_AUC, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_color_manual(values = colors_) +  # Specify the custom color palette+
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of point sizes as needed
  ylim(0.38, 1.02) +
  theme_minimal()+
  labs(x = "-log10 of p-values HLA-CMV association ", y = "Specific AUC", title = "HLA-I", color = "Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))


p2 = ggplot(df2, aes(x = -log10(HLA_CMV_asso_pval) , y = specific_AUC, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_hline(yintercept = 0.8,linetype = "dashed", color = "red") +
  scale_color_manual(values = colors_) +  # Specify the custom color palette+
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of point sizes as needed
  ylim(0.38, 1.02) +
  theme_minimal()+
  labs(x = "-log10 of p-values HLA-CMV association ", y = "Specific AUC", title = "HLA-II", color="Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))

legend1 <- get_legend(p1)
legend2 <- get_legend(p2)

p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

pdf("figure4_panels/Fig4A.pdf", width = 4.8, height = 3.2)
grid.arrange(
  arrangeGrob(
    p1, p2,
    legend1,legend2,
    nrow=2,
    heights = c(56, 4)  # Adjust heights to allocate space for the legend
    #top = textGrob("HLA-I HLA-II Plots", gp = gpar(fontsize = 14, fontface = "bold"))
  ),
  bottom = textGrob("-log10 of p-values HLA-CMV association", gp = gpar(fontsize = 12)),
  left = textGrob("HLA-specific AUC", rot = 90, gp = gpar(fontsize = 12))
)
dev.off()

#######################
### regression plot ###
#######################

# remove an outlier
df2 = df2 %>% filter(specific_training_size < 200) 
HLAI_fit <- rq(specific_AUC ~ specific_training_size, data = df1, tau = 0.5)
HLAII_fit = rq(specific_AUC ~ specific_training_size, data = df2, tau = 0.5)
df1$regression_pred = predict(HLAI_fit)
df2$regression_pred = predict(HLAII_fit)


reg1 = ggplot(df1, aes(x = specific_training_size, y = specific_AUC, color = HLA_group)) +
  geom_point(size=2.5, alpha = 0.6, show.legend = T) +
  geom_line(aes(y = regression_pred), color = "blue") + # Fitted quantile regression line
  geom_hline(yintercept = 0.80,linetype = "dashed", color = "red") +
  scale_color_manual(name = "Group", values = colors_) +
  ylim(0.38, 1.02) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))


reg2 = ggplot(df2, aes(x = specific_training_size, y = specific_AUC, color = HLA_group)) +
  geom_point(size=2.5, alpha = 0.6, show.legend = T) +
  geom_line(aes(y = regression_pred), color = "blue") + # Fitted quantile regression line
  geom_hline(yintercept = 0.8,linetype = "dashed", color = "red") +
  scale_color_manual(name = "Group", values = colors_) +
  ylim(0.38, 1.02) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))


legend1 <- get_legend(reg1)
legend2 <- get_legend(reg2)

reg1 <- reg1 + ggtitle("HLA-I") +  theme(legend.position = "none", plot.title = element_text(size=14, face="bold", hjust = 0.5)) 
reg2 <- reg2 +  ggtitle("HLA-II") + theme(legend.position = "none", plot.title = element_text(size=14, face="bold", hjust = 0.5)) 

pdf("figure4_panels/Fig4B.pdf", width = 4.8, height = 3.2)
grid.arrange(
  arrangeGrob(
    reg1, reg2,
    legend1,legend2,
    nrow=2,
    heights = c(56, 4) # Adjust heights to allocate space for the legend
  ),
  bottom = textGrob("Training sample size", gp = gpar(fontsize = 12)),
  left = textGrob("HLA-specific AUC", rot = 90, gp = gpar(fontsize = 12))
)
dev.off()


sessionInfo()
q(save = "no")


