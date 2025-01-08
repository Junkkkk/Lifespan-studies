#!/bin/bash

### All
## Make GRM
gcta64 \
--bfile ukb_wes_merged_393833_prune_maf_0.01 \
--autosome \ 
--make-grm \
--out ukb_wes_merged_393833_prune_maf_0.01 \
--threads 40

## Make sparse matrix
gcta64 \
--grmukb_wes_merged_393833_prune_maf_0.01 \
--make-bK-sparse 0.05 \
--out ukb_wes_merged_393833_prune_maf_0.01_sparse \
--threads 40

## Test mixed model
gcta64 \
--grm-sparse ukb_wes_merged_393833_prune_maf_0.01_sparse \
--fastGWA-mlm \
--pheno pheno_martingale_residuals.csv \
--bfile ukb_wes_merged_393833_prune_maf_0.01 \
--save-fastGWA-mlm-residual \
--out adjusted_martingale_residuals \
--threads 40

### Male
## Make GRM
gcta64 \
--bfile ukb_wes_merged_393833_prune_maf_0.01_male \
--autosome \ 
--make-grm \
--out ukb_wes_merged_393833_prune_maf_0.01_male \
--threads 40

## Make sparse matrix
gcta64 \
--grmukb_wes_merged_393833_prune_maf_0.01_male \
--make-bK-sparse 0.05 \
--out ukb_wes_merged_393833_prune_maf_0.01_sparse_male \
--threads 40

## Test mixed model
gcta64 \
--grm-sparse ukb_wes_merged_393833_prune_maf_0.01_sparse_male \
--fastGWA-mlm \
--pheno pheno_martingale_residuals_male.csv \
--bfile ukb_wes_merged_393833_prune_maf_0.01_male \
--save-fastGWA-mlm-residual \
--out adjusted_martingale_residuals_male \
--threads 40

### Female
## Make GRM
gcta64 \
--bfile ukb_wes_merged_393833_prune_maf_0.01_female \
--autosome \ 
--make-grm \
--out ukb_wes_merged_393833_prune_maf_0.01_female \
--threads 40

## Make sparse matrix
gcta64 \
--grmukb_wes_merged_393833_prune_maf_0.01_female \
--make-bK-sparse 0.05 \
--out ukb_wes_merged_393833_prune_maf_0.01_sparse_female \
--threads 40

## Test mixed model
gcta64 \
--grm-sparse ukb_wes_merged_393833_prune_maf_0.01_sparse_female \
--fastGWA-mlm \
--pheno pheno_martingale_residuals_male.csv \
--bfile ukb_wes_merged_393833_prune_maf_0.01_female \
--save-fastGWA-mlm-residual \
--out adjusted_martingale_residuals_female \
--threads 40