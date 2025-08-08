#

df1 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA I")
df2 = read_excel("../../codes/result_tables/Table_S1.xlsx", sheet = "HLA II")

df1$HLA[which(df1$CMV_asso_pvalue<0.05)]
df2$HLA[which(df2$CMV_asso_pvalue<0.05)]