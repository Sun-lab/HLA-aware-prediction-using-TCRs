library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)

setwd("/home/hyo/mystuff/git_test/conditional_TCR_prediction/plot/figure3/")

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

p1 = ggplot(df1, aes(x = HLA_specific_AUC, y = Combined_AUC, size = Training_Size, color = HLA_group)) +
  geom_point(alpha = 0.6, show.legend = T) +
  scale_size(name = "Size", breaks = c(50,100,150), labels=c("50","100","150") )+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") + 
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") + 
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") + 
  scale_color_manual(name = "   Group", values = colors_) +  # Specify the custom color palette+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = "none", 
         size = guide_legend(nrow = 1))


p2 = ggplot(df1, aes(x = HLA_ignorant_AUC, y = Combined_AUC, size = Training_Size, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, show.legend = T) +
  scale_size(name = "Size", breaks = c(50,100,150), labels=c("50","100","150") )+
  geom_text_repel(data=subset(df1, abs(Combined_AUC - HLA_ignorant_AUC)>0.05 ), size = 1.5, vjust = 1, show.legend = FALSE) +  # Use ggrepel for better label placement
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") + 
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") + 
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") + 
  scale_color_manual(name = "   Group", values = colors_) + 
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = "none", 
         size = guide_legend(nrow = 1))

p3 = ggplot(df2, aes(x = HLA_specific_AUC, y = Combined_AUC, size = Training_Size, color = HLA_group)) +
  geom_point(alpha = 0.6, show.legend = T) +
  scale_size(name = "Size", breaks = c(50,100,150), labels=c("50","100","150") )+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") +  
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") +  
  scale_color_manual(name = "   Group", values = colors_) +  # Specify the custom color palette+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = "none", 
         size = guide_legend(nrow = 1))

p4 = ggplot(df2, aes(x = HLA_ignorant_AUC, y = Combined_AUC, size = Training_Size, color = HLA_group, label = HLA_name)) +
  geom_point(alpha = 0.6, show.legend = T) +
  scale_size(name = "Size", breaks = c(50,100,150), labels=c("50","100","150") )+
  geom_text_repel(data=subset(df2, abs(Combined_AUC - HLA_ignorant_AUC)>0.03), size = 1, vjust = -1, show.legend = FALSE) +  # Use ggrepel for better label placement
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = baseline,linetype = "dashed", color = "red") +  
  geom_hline(yintercept = baseline,linetype = "dashed", color = "red") +  
  scale_color_manual(name = "   Group", values = colors_) +  
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(0.01, "mm") ) +
  guides(color = "none", 
         size = guide_legend(nrow = 1))


grid.arrange(
  arrangeGrob(
    p1, p2,
    p3,p4,
    nrow=1
  ),
  bottom = textGrob("Training Size", gp = gpar(fontsize = 12)),
  left = textGrob("HLA-specific AUC", rot = 90, gp = gpar(fontsize = 12))
)
