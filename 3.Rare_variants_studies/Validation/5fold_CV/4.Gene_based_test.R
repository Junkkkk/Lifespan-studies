library(data.table)
library(SKAT)


for (i in 1:5) {
  
  pheno_txt <- paste0('pheno_393833_all_fold', i, '.txt')
  completeID_txt <- paste0('completeDataID_393833_all_fold_', i, '.txt')
  pheno_covar_file <- paste0('pheno_covar_393833_caucasian_survial_death_residuals_fold', i, '.csv')
  
  pheno_covar <- read.csv(pheno_covar_file)
  
  ### extract bfile ###
  anno <- fread(anno_csv_file)
  
  if (nrow(anno) == 0){
  return(NULL)
  }
  
  # 1, 2, 3, ... , 22
  chr_num <- gsub(".*/chr(\\d+)/.*", "\\1", anno_csv_file)
  # APOE, BRCA1, ...
  gene <- gsub(".csv","",strsplit(basename(anno_csv_file), "_")[[1]][1])
  
  if (file.exists(file.path(skato_result_path, 'all_csv', paste0('chr', chr_num, '_', gene,'.csv')))) {
  return(NULL)
  }
  
  ifelse(!dir.exists(gene_path), dir.create(gene_path), FALSE)
  ifelse(!dir.exists(file.path(gene_path, paste0('chr',chr_num, '_',gene))), dir.create(file.path(gene_path, paste0('chr',chr_num, '_',gene))), FALSE)
  
  write.table(anno$SNP, file = file.path(gene_path, paste0('chr',chr_num,'_',gene), 'extract.csv'), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  system(paste('plink2 --bfile', paste0(bfile_path, '/ukb_wes_chr', chr_num),
               '--allow-no-sex --extract', file.path(gene_path, paste0('chr',chr_num,'_',gene), 'extract.csv'),
               '--keep', file.path(pheno_txt_path, completeID_txt),
               '--pheno', file.path(pheno_txt_path, pheno_txt),
               '--make-bed --out', file.path(gene_path, paste0('chr',chr_num,'_',gene), paste0('ukb_wes_chr', chr_num, '_', gene))
               )
         )
  
  ### gene-based test ###
  ifelse(!dir.exists(skato_result_path), dir.create(skato_result_path), FALSE)
  ifelse(!dir.exists(file.path(skato_result_path,'all_csv')), dir.create(file.path(skato_result_path,'all_csv')), FALSE)
  
  gene_folder <- file.path (gene_path, paste0('chr', chr_num, '_', gene))
  
  File.Bed <- paste0(gene_folder, '/ukb_wes_chr',chr_num,'_',gene, '.bed')
  File.Bim <- paste0(gene_folder, '/ukb_wes_chr',chr_num,'_',gene, '.bim')
  File.Fam <- paste0(gene_folder, '/ukb_wes_chr',chr_num,'_',gene, '.fam')
  
  fam <- read.table(File.Fam, sep='\t', header=F)
  fam <- fam[,1:2]
  colnames(fam) <- c("FID","IID")
  fam$order <- 1:nrow(fam)
  
  dfm <- merge(fam, pheno_covar, by=c('FID','IID'))
  dfm <- dfm[order(dfm$order), ]
  
  # Definition of variants
  snps_lof <- anno[anno$IMPACT == 'HIGH', 'SNP']
  snps_alpha_missense_70 <- anno[(grepl('likely_pathogenic', anno$am_class) & anno$am_pathogenicity>=0.7), 'SNP']
  snps_missense_REVEL_75 <- anno[(grepl('missense_variant', anno$Consequence) & anno$REVEL>=0.75), 'SNP']
  
  SNPs <- list(snps_lof, snps_alpha_missense_70, snps_missense_REVEL_75)
  names <- c('lof', 'alpha_missense_70', 'missense_REVEL_75')
  results <- data.frame('gene', 'SNPs_category', 'SKATO_p', 'N_snps', 'N_snps_test','method')
  
  if (sum(sapply(SNPs, nrow))==0){
    return(NULL)
  }
  
  for (i in 1:length(SNPs)){
    if (nrow(SNPs[[i]]) > 0){
      df1 <- data.frame(Gene=replicate(nrow(SNPs[[i]]), gene), SNP=SNPs[[i]])
  
      # write list of considered variants to file
      File.SetID <- paste0(gene_folder, '/', names[i], '.SetID')
      write.table(df1, File.SetID, col.names=F, row.names=F, quote=F)
  
      # generate SSD file for SKAT
      File.SSD <- paste0(gene_folder, '/', names[i], '.SSD')
      File.Info <- paste0(gene_folder, '/', names[i], '.Info')
      Generate_SSD_SetID(File.Bed,File.Bim,File.Fam,File.SetID,File.SSD,File.Info)
  
      # open SSD file
      SSD.Info <- Open_SSD(File.SSD, File.Info)
      Set.Index <- 1 # only one gene set in this SSD file
      Z <- Get_Genotypes_SSD(SSD.Info, Set.Index, is_ID = F)
  
      # Burden
      prelim <- SKAT_Null_Model(dfm$adjusted_martingale_residuals_fold ~ 1, out_type="C")
      fit <- SKAT(Z, prelim, method='Burden', kernel="linear.weighted", weights.beta=c(1,25), is_dosage = FALSE, estimate_MAF=1)
      results <- rbind(results, c(gene, names[i], fit$p.value, fit$param$n.marker, fit$param$n.marker.test, 'Burden'))
  
      # SKAT-O
      prelim <- SKAT_Null_Model(dfm$adjusted_martingale_residuals_fold ~ 1, out_type="C")
      fit <- SKAT(Z, prelim, method='SKATO', kernel="linear.weighted", weights.beta=c(1,25), is_dosage = FALSE, estimate_MAF=1)
      results <- rbind(results, c(gene, names[i], fit$p.value, fit$param$n.marker, fit$param$n.marker.test, 'SKATO'))
  
      # close SSD file
      Close_SSD()
    }
  }
  t <- data.frame(results)[-1, ]
  colnames(t) <- c('gene', 'SNPs_category', 'SKATO_p', 'N_snps', 'N_snps_test','method')
  
  if (length(t)>0){
  write.csv(t, 'result.csv',row.names=F,quote=F)
  }
}