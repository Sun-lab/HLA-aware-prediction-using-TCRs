library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)
library(quantreg)



setwd("/home/hyo/mystuff/git_test/conditional_TCR_prediction/plot/figure4/")

get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


df1 = read.csv("../HLA_I_ready_data.csv")
df2 = read.csv("../HLA_II_ready_data.csv")
baseline = 0.9163461538461538
colors_ <- c("#1f77b4", "#ff7f0e", "#2ca02c")

########################
# Specific-AUC v.s. pval
########################

p1 = ggplot(df1, aes(x = -log10(HLA_CMV_asso_pval) , y = HLA_specific_AUC, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_hline(yintercept = 0.796, linetype = "dashed", color = "red") +
  scale_color_manual(values = colors_) +  # Specify the custom color palette+
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of point sizes as needed
  theme_minimal()+
  labs(x = "-Log10 of p-values HLA-CMV association ", y = "Specific AUC", title = "HLA-I", color = "Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))


p2 = ggplot(df2, aes(x = -log10(HLA_CMV_asso_pval) , y = HLA_specific_AUC, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_hline(yintercept = 0.8,linetype = "dashed", color = "red") +
  scale_color_manual(values = colors_) +  # Specify the custom color palette+
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of point sizes as needed
  theme_minimal()+
  labs(x = "-Log10 of p-values HLA-CMV association ", y = "Specific AUC", title = "HLA-II", color="Group") +
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

grid.arrange(
  arrangeGrob(
    p1, p2,
    legend1,legend2,
    nrow=2,
    heights = c(56, 4)  # Adjust heights to allocate space for the legend
    #top = textGrob("HLA-I HLA-II Plots", gp = gpar(fontsize = 14, fontface = "bold"))
  ),
  bottom = textGrob("-Log10 of p-values HLA-CMV association", gp = gpar(fontsize = 12)),
  left = textGrob("HLA-specific AUC", rot = 90, gp = gpar(fontsize = 12))
)


#######################
### regression plot ###
#######################

# remove an outlier
df2 = df2 %>% filter(Training_Size < 200) 
HLAI_fit <- rq(HLA_specific_AUC ~ Training_Size, data = df1, tau = 0.5)
HLAII_fit = rq(HLA_specific_AUC ~ Training_Size, data = df2, tau = 0.5)
df1$regression_pred = predict(HLAI_fit)
df2$regression_pred = predict(HLAII_fit)




reg1 = ggplot(df1, aes(x = Training_Size, y = HLA_specific_AUC, color = HLA_group)) +
  geom_point(size=2.5, alpha = 0.6, show.legend = T) +
  geom_line(aes(y = regression_pred), color = "blue") + # Fitted quantile regression line
  geom_hline(yintercept = 0.79,linetype = "dashed", color = "red") +
  scale_color_manual(name = "Group", values = colors_) +
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))




reg2 = ggplot(df2, aes(x = Training_Size, y = HLA_specific_AUC, color = HLA_group)) +
  geom_point(size=2.5, alpha = 0.6, show.legend = T) +
  geom_line(aes(y = regression_pred), color = "blue") + # Fitted quantile regression line
  geom_hline(yintercept = 0.8,linetype = "dashed", color = "red") +
  scale_color_manual(name = "Group", values = colors_)
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


grid.arrange(
  arrangeGrob(
    reg1, reg2,
    legend1,legend2,
    nrow=2,
    heights = c(56, 4) # Adjust heights to allocate space for the legend
  ),
  bottom = textGrob("Training Size", gp = gpar(fontsize = 12)),
  left = textGrob("HLA-specific AUC", rot = 90, gp = gpar(fontsize = 12))
)




