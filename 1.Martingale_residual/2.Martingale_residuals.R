library(survival)

df <- read.csv('pheno.csv')

## All ##

#PC1-40
pc_vars <- paste("PC", 1:40, sep = "")
#geographic covariates as dummy variables
c_vars <- paste("c", 1:15, sep = "")

formula <- as.formula(paste("Surv(Age_1st_visit, last_known_age, Death) ~ Sex +", paste(pc_vars, collapse = " + "), "+", paste(c_vars, collapse = " + ")))

# Fit the Cox model
fit <- coxph(formula, data = df)

martingale_residuals_all <- residuals(fit, type = 'martingale')
merged_data <- cbind(df, residuals)

write.csv(merged_data, 'pheno_martingale_residuals.csv', quote=F, row.names=F)

## Male ##

#PC1-40
df <- read.csv('pheno.csv')
df <- df[df$Sex==1,]

pc_vars <- paste("PC", 1:40, sep = "")
#geographic covariates as dummy variables
c_vars <- paste("c", 1:15, sep = "")

formula <- as.formula(paste("Surv(Age_1st_visit, last_known_age, Death) ~", paste(pc_vars, collapse = " + "), "+", paste(c_vars, collapse = " + ")))

# Fit the Cox model
fit <- coxph(formula, data = df)

martingale_residuals_male <- residuals(fit, type = 'martingale')
merged_data <- cbind(df, residuals)

write.csv(merged_data, 'pheno_martingale_residuals_male.csv', quote=F, row.names=F)

## Female ##

#PC1-40
df <- read.csv('pheno.csv')
df <- df[df$Sex==0,]

pc_vars <- paste("PC", 1:40, sep = "")
#geographic covariates as dummy variables
c_vars <- paste("c", 1:15, sep = "")

formula <- as.formula(paste("Surv(Age_1st_visit, last_known_age, Death) ~", paste(pc_vars, collapse = " + "), "+", paste(c_vars, collapse = " + ")))

# Fit the Cox model
fit <- coxph(formula, data = df)

martingale_residuals_female <- residuals(fit, type = 'martingale')
merged_data <- cbind(df, residuals)

write.csv(merged_data, 'pheno_martingale_residuals_female.csv', quote=F, row.names=F)