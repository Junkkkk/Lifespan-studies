library(survival)

df <- read.csv("pheno_covar_73281_independent_samples.csv")

# Define variables for principal components (PCs) and geographic covariates
pc_vars <- paste("PC", 1:40, sep = "")  # Principal components (PC1 to PC40)
c_vars <- paste("c", 1:15, sep = "")   # Geographic covariates (c1 to c15)

# Construct the formula for the Cox proportional hazards model including ancestry dummy variables
formula <- as.formula(paste(
  "Surv(Age_1st_visit, last_known_age, Death) ~ Sex + Black + Mixed + Other + White +",
  paste(pc_vars, collapse = " + "), "+",
  paste(c_vars, collapse = " + ")
))

# Fit the Cox proportional hazards model
fit <- coxph(formula, data = df)

# Calculate martingale residuals for each individual
martingale_residuals <- residuals(fit, type = "martingale")

merged_data <- cbind(df, martingale_residuals)

write.csv(
  merged_data, 
  paste0("pheno_martingale_residuals_independent_samples.csv"), 
  quote = FALSE, 
  row.names = FALSE
)