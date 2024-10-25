import os


bolt = os.path.join(DIR_GWAS_QC, 'BOLT-LMM_v2.3.4', 'bolt')
ldScoresFile = os.path.join(DIR_GWAS_QC, 'BOLT-LMM_v2.3.4', 'tables', 'LDSCORE.1000G_EUR.tab.gz')
geneticMap = os.path.join(DIR_GWAS_QC, 'BOLT-LMM_v2.3.4', 'tables', 'genetic_map_hg19_withX.txt.gz')

phenCol = 'martingale_residuals'

cmd = ' '.join([bolt,
                '--bed=%s' % DIR_UKB_CALL+'/ukb_cal_chr{1:22}_v2.bed',
                '--bim=%s' % DIR_UKB_CALL+'/ukb_cal_chr{1:22}_v2.bim',
                '--fam=%s' % fam0,
                '--bgenFile=%s' % DIR_UKB_IMP+'/ukb_imp_chr{1:22}_v3.bgen',
                '--bgenMinMAF=1e-3',
                '--bgenMinINFO=0.3',
                '--remove=%s' % remove1,
                '--phenoFile=%s' % pheno,
                '--phenoCol=%s' % phenCol,
                '--LDscoresFile=%s' % ldScoresFile,
                '--geneticMapFile=%s' % geneticMap,
                '--lmmForceNonInf',
                '--numThreads=30',
                '--sampleFile=%s' % sample,  
                '--statsFile=%s' % output,
                '--statsFileBgenSnps=%s' % output_gen,
                '--verboseStats'])

print (cmd)
os.system(cmd)