library(data.table)
library(survival)
library(tidyverse)

pheno_covar <- fread('pheno_covar.csv')
res <- fread('all_results.csv')
res <- res[res$P < 7.4e-07]

genes_variant <- data.frame(gene = res$gene, variant_type = res$variant_type)

results <- data.frame(chr = character(), gene = character(), variant_type = character(), HR = numeric(), P = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(genes_variant)) {
  
  gene <- genes_variant$gene[i]
  variant_type <- genes_variant$variant_type[i]
  chm <- paste0('chr', res[res$gene == gene, 'chr'][1])
  
  SNP_ID <- fread(paste0(chm, '_', gene, '/', variant_type, '.SetID'), header = FALSE)
  fwrite(SNP_ID[, 2], file = paste0(chm, '_', gene, '/', variant_type, '_SNP_ID.txt'), col.names = FALSE)
  
  system(paste0('plink --bfile ', chm, '_', gene, '/ukb_wes_', chm, '_', gene,
                ' --extract ', chm, '_', gene, '/', variant_type, '_SNP_ID.txt',
                ' --recode A --out ukb_wes_', chm, '_', gene, '_', variant_type))
  
  
  
  raw_file <- fread(paste0(chm, '_', gene, '/ukb_wes_', chm, '_', gene, '_', variant_type, '.raw'), sep = ' ')
  
  # Define burden
  row_sums <- rowSums(raw_file[, 7:ncol(raw_file)])
  raw_file[[paste0(variant_type, '_carrier')]] <- ifelse(row_sums != 0, 1, 0)
  
  df <- merge(pheno_covar, raw_file, by = c('FID', 'IID'))
  
  # Add covariate
  pc_columns <- paste0("PC", 1:40)
  c_columns <- paste0("c", 1:15)
  
  formula <- as.formula(paste("Surv(last_known_age, Death) ~ Sex + ", paste(pc_columns, collapse = " + "), " + ",
                              paste(c_columns, collapse = " + "), " + ", variant_type, "_carrier + cluster(FID)"))
  
  # Run Cox regression model
  cox_model <- coxph(formula, data = df)
  summary_model <- summary(cox_model)
  
  p_value <- summary_model$coefficients[paste0(variant_type, "_carrier"), "Pr(>|z|)"]
  hr <- exp(coef(cox_model)[paste0(variant_type, "_carrier")])
  rounded_hr <- round(hr, 1)
  p_value <- format(p_value, scientific = TRUE, digits = 1)
  
  info <- data.frame(chr = chm, gene = gene, variant_type = variant_type, HR = rounded_hr, P = p_value)
  results <- rbind(results, info)
}

print(results)