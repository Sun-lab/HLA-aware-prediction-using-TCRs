#

df1 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA I")
df2 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA II")

df1$diff = df1$combined_AUC - df1$agnostic_AUC
df2$diff = df2$combined_AUC - df2$agnostic_AUC

summary(df1$diff)
summary(df2$diff)

wilcox.test(df1$combined_AUC, df1$agnostic_AUC, paired = TRUE, alternative = "greater")
wilcox.test(df2$combined_AUC, df2$agnostic_AUC, paired = TRUE, alternative = "greater")

wilcox.test(df1$diff, mu = 0, alternative = "greater")
wilcox.test(df2$diff, mu = 0, alternative = "greater")

df1 = df1[which(df1$run_10splits==FALSE),]
df2 = df2[which(df2$run_10splits==FALSE),]

summary(df1$diff)
summary(df2$diff)

wilcox.test(df1$combined_AUC, df1$agnostic_AUC, paired = TRUE, alternative = "greater")
wilcox.test(df2$combined_AUC, df2$agnostic_AUC, paired = TRUE, alternative = "greater")

wilcox.test(df1$diff, mu = 0, alternative = "greater")
wilcox.test(df2$diff, mu = 0, alternative = "greater")
