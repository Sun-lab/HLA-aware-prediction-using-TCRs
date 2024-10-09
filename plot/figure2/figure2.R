library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)



setwd("/home/hyo/mystuff/git_test/conditional_TCR_prediction/plot/figure2/")
df1 = read.csv("../HLA_I_ready_data.csv")
df2 = read.csv("../HLA_II_ready_data.csv")

baseline = 0.9163461538461538
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c")

get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## (A)

f1_main = ggplot(df1, aes(x = HLA_ignorant_AUC, y = HLA_specific_AUC, size = Training_Size, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") +  
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") +  
  scale_color_manual(labels =c("A/DP","B/DQ","C/DR") , values = my_colors) +  
  scale_size(range = c(4,8), breaks = c(50,100,150))+
  labs(x = "all AUC", y = "specific AUC", size = "train size", title = "HLA-I") +
  theme_minimal()  +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = guide_legend(nrow = 1), 
         size = "none")



f1_main_2 = ggplot(df2, aes(x = HLA_ignorant_AUC, y = HLA_specific_AUC, size = Training_Size, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, show.legend = T) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") +  
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") +  
  scale_color_manual(values = my_colors) +  
  scale_size(range = c(4,8), breaks = c(50,100,150))+
  labs(x = "all AUC", y = "specific AUC", color = "group", size="size", title = "HLA-II") +
  theme_minimal()   +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(color = "none",
         size = guide_legend(nrow = 1))


legend1 <- get_legend(f1_main)
legend2 <- get_legend(f1_main_2)

f1_main <- f1_main + theme(legend.position = "none")
f1_main_2 <- f1_main_2 + theme(legend.position = "none")

grid.arrange(
  arrangeGrob(
    f1_main, f1_main_2,
    legend1,legend2,
    nrow=2,
    heights = c(35, 4)  # Adjust heights to allocate space for the legend
  ),
  bottom = textGrob("All model AUC", gp = gpar(fontsize = 10)),
  left = textGrob("Specific model AUC", rot = 90, gp = gpar(fontsize = 12))
)


# (B)
par(mfrow=c(1,2))
hist(df1$Significant_TCR_prop, main="HLA-I", xlab = "")
hist(df2$Significant_TCR_prop, main="HLA-II", xlab ="")


# (C) 

i3_p1 = ggplot(df1, aes(x = log(Training_Size), y = log(Significant_TCR_prop*280), color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, show.legend = T, size=2) +
  geom_text_repel(data=subset(df1, Significant_TCR_prop*280 >= 60 ), size = 3, vjust = -1, show.legend = FALSE) +  
  scale_color_manual(values = my_colors) +  
  scale_size_continuous(range = c(3, 10)) +  
  labs(x = "all AUC", y = "specific AUC", size = "train size", title = "HLA-I") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))

i3_p2 = ggplot(df2, aes(x = log(Training_Size), y = log(Significant_TCR_prop*280), color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, show.legend = T, size=2) +
  geom_text_repel(data=subset(df2, Significant_TCR_prop*280 > 60 ),size = 3, vjust = -1, show.legend = FALSE) +  
  scale_color_manual(values = my_colors) +  
  scale_size_continuous(range = c(3, 10)) +  
  labs(x = "all AUC", y = "specific AUC", size = "train size", title = "HLA-II") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1))

legend1 <- get_legend(i3_p1)
legend2 <- get_legend(i3_p2)

i3_p1 <- i3_p1 + theme(legend.position = "none")
i3_p2 <- i3_p2 + theme(legend.position = "none")

grid.arrange(
  arrangeGrob(
    i3_p1, i3_p2,
    legend1,legend2,
    nrow=2,
    heights = c(56, 4)  # Adjust heights to allocate space for the legend
  ),
  bottom = textGrob("Log10 of Traning Sample Size", gp = gpar(fontsize = 12)),
  left = textGrob("Log10 of Number of Significant TCRs", rot = 90, gp = gpar(fontsize = 12))
)

