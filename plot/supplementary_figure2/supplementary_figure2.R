library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)
library(readxl)

df1 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA I")
df2 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA II")

df1$size_category = rep("<=70", nrow(df1))
df1$size_category[which(df1$HLA_frequency>70)] = ">70"

df2$size_category = rep("<=70", nrow(df2))
df2$size_category[which(df2$HLA_frequency>70)] = ">70"

baseline = 0.9163461538461538
colors_ <- c("#1f77b4", "#ff7f0e", "#2ca02c")

p1 = ggplot(df1, aes(x = specific_AUC, y = knn_AUC, size = size_category, color = HLA_group)) +
  geom_point(alpha = 0.6, show.legend = T) +
  scale_size_manual(
    name = "Freq",
    values = c("<=70" = 2, ">70" = 3),  # adjust sizes as needed
    labels = c("<=70", ">70")
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") + 
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") + 
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") + 
  scale_color_manual(labels =c("A/DP","B/DQ","C/DR") , values = colors_, name = "group")+  # Specify the custom color palette+
  xlab("HLA-specific AUC") + 
  ylab("KNN AUC") + 
  xlim(0.38, 1.02) + 
  ylim(0.36, 1.02) +
  theme_minimal()+
  ggtitle("HLA-I") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),  # trim outer margin
        legend.spacing.y = unit(0, "mm")) +
  guides(color = guide_legend(nrow = 1), 
         size = "none")

p2 = ggplot(df2, aes(x = specific_AUC, y = knn_AUC, size = size_category, color = HLA_group)) +
  geom_point(alpha = 0.6, show.legend = T) +
  scale_size_manual(
    name = "Freq",
    values = c("<=70" = 2, ">70" = 3),  # adjust sizes as needed
    labels = c("<=70", ">70")
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") +  
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") +  
  scale_color_manual(labels =c("A/DP","B/DQ","C/DR") , values = colors_, name = "group") +  # Specify the custom color palette+
  xlab("HLA-specific AUC") + 
  ylab("KNN AUC") + 
  xlim(0.38, 1.02) + 
  ylim(0.36, 1.02) +
  theme_minimal()+
  ggtitle("HLA-II") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),  # trim outer margin
        legend.spacing.y = unit(0, "mm") ) +
  guides(color = "none",
         size = guide_legend(nrow = 1))

pdf("supplementary_figure2_panels/Supplementary_Fig2.pdf", width = 5.2, height = 3.2)
grid.arrange(
  arrangeGrob(
    p1, p2,
    nrow=1
  )
)
dev.off()


sessionInfo()
q(save = "no")


