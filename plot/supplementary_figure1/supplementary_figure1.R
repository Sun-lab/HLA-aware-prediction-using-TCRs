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



df_hla_names = read.csv("../../codes/intermediate_files/complete_HLA_rownames.csv",
                        header=TRUE)
hla_names = df_hla_names$HLA_name

#######################
### regression plot ###
#######################

# keep the outlier
# df2 = df2 %>% filter(specific_training_size < 200) 
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
  ylim(0.38, 1.7) +
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

pdf("supplementary_figure1_panels/Supplementary_Fig1.pdf", width = 4.8, height = 3.2)
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


