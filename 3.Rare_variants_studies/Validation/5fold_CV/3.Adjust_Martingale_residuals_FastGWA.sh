#!/bin/bash


# Loop through folds 1 to 5
for i in {1..5}; do

  ## Step 1: PLINK
  plink \
    --bfile ukb_wes_merged_393833_prune_maf_0.01 \
    --keep completeDataID_393833_all_fold${i}.txt \
    --out ukb_wes_merged_393833_prune_maf_0.01_fold${i}

  ## Step 2: Make GRM
  gcta64 \
    --bfile ukb_wes_merged_393833_prune_maf_0.01_fold${i} \
    --autosome \
    --make-grm \
    --out ukb_wes_merged_393833_prune_maf_0.01_fold${i} \
    --threads 40

  ## Step 3: Make sparse matrix
  gcta64 \
    --grm ukb_wes_merged_393833_prune_maf_0.01_fold${i} \
    --make-bK-sparse 0.05 \
    --out ukb_wes_merged_393833_prune_maf_0.01_fold${i}_sparse \
    --threads 40

  ## Step 4: Test mixed model
  gcta64 \
    --grm-sparse ukb_wes_merged_393833_prune_maf_0.01_fold${i}_sparse \
    --fastGWA-mlm \
    --pheno pheno_martingale_residuals_fold${i}.csv \
    --bfile ukb_wes_merged_393833_prune_maf_0.01_fold${i} \
    --save-fastGWA-mlm-residual \
    --out adjusted_martingale_residuals_fold${i} \
    --threads 40

done