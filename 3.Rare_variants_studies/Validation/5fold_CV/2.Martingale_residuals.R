library(survival)


for (i in 1:5) {
  df <- read.csv(paste0("pheno_covar_fold", i, ".csv"))
  
  # Define variables for PCs and geographic covariates
  pc_vars <- paste("PC", 1:40, sep = "")
  c_vars <- paste("c", 1:15, sep = "")
  
  # Construct the formula for the Cox model
  formula <- as.formula(paste("Surv(Age_1st_visit, last_known_age, Death) ~ Sex +",
                              paste(pc_vars, collapse = " + "), "+",
                              paste(c_vars, collapse = " + ")))
  
  # Fit the Cox proportional hazards model
  fit <- coxph(formula, data = df)
  
  # Calculate martingale residuals
  martingale_residuals_fold <- residuals(fit, type = "martingale")
  
  # Merge the residuals with the original data
  merged_data <- cbind(df, martingale_residuals_fold)
  
  # Write the merged data to a new CSV file
  write.csv(merged_data, paste0("pheno_martingale_residuals_fold", i, ".csv"), 
            quote = FALSE, row.names = FALSE)
}