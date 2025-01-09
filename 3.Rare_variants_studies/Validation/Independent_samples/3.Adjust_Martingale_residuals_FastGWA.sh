#!/bin/bash

  ## Step 1: PLINK
  plink \
    --bfile All_individuals_kept_quality_score_0.3_geno_0.05_maf_0.01_hwe_1e-6_pruned \
    --keep completeDataID_73281.txt \
    --out ukb_wes_merged_73281_prune_maf_0.01

  ## Step 2: Make GRM
  gcta64 \
    --bfile ukb_wes_merged_73281_prune_maf_0.01 \
    --autosome \
    --make-grm \
    --out ukb_wes_merged_73281_prune_maf_0.01 \
    --threads 40

  ## Step 3: Make sparse matrix
  gcta64 \
    --grm ukb_wes_merged_73281_prune_maf_0.01 \
    --make-bK-sparse 0.05 \
    --out ukb_wes_merged_73281_prune_maf_0.01_sparse \
    --threads 40

  ## Step 4: Test mixed model
  gcta64 \
    --grm-sparse ukb_wes_merged_73281_prune_maf_0.01_sparse \
    --fastGWA-mlm \
    --pheno pheno_martingale_residuals_independent_samples.csv \
    --bfile ukb_wes_merged_73281_prune_maf_0.01 \
    --save-fastGWA-mlm-residual \
    --out adjusted_martingale_residuals_independent_samples \
    --threads 40