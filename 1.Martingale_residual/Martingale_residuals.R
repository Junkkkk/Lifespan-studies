library(survival)

df <- read.csv('pheno.csv')

#PC1-40
pc_vars <- paste("PC", 1:40, sep = "")
#geographic covariates as dummy variables
c_vars <- paste("c", 1:15, sep = "")

formula <- as.formula(paste("Surv(Age_1st_visit, last_known_age, Death) ~ Sex +", paste(pc_vars, collapse = " + "), "+", paste(c_vars, collapse = " + ")))

# Fit the Cox model
fit <- coxph(formula, data = df)

residuals <- residuals(fit, type = 'martingale')
merged_data <- cbind(df, residuals)

write.csv(merged_data, 'pheno_martingale_residuals.csv', quote=F, row.names=F)