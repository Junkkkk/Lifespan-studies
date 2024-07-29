library(survival)

df <- read.csv('pheno.csv')

fit <- coxph(Surv(last_known_Age, Death) ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5, data = df)

residuals <- residuals(fit, type = 'martingale')
merged_data <- cbind(df, residuals)

write.csv(merged_data, 'pheno_martingale_residuals.csv', quote=F, row.names=F)