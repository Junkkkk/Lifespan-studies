import os

### group : 'all / male / female'
group='all'

fam0 = f'ukb_cal_chr1_v2_{group}.fam'
pheno = f'pheno_martingale_residuals_{subgroupjects}.csv'
phenCol = f'martingale_residuals_{group}'

output = os.path.join('stats_'+phenCol)
output_gen = os.path.join('stats_bgen_'+phenCol+f'_{group}')

bolt = os.path.join(DIR_GWAS_QC, 'BOLT-LMM_v2.3.4', 'bolt')
ldScoresFile = os.path.join(DIR_GWAS_QC, 'BOLT-LMM_v2.3.4', 'tables', 'LDSCORE.1000G_EUR.tab.gz')
geneticMap = os.path.join(DIR_GWAS_QC, 'BOLT-LMM_v2.3.4', 'tables', 'genetic_map_hg19_withX.txt.gz')
sample = 'ukb45420_imp_chrX_v3_s486669.sample'

cmd = ' '.join([bolt,
                '--bed=%s' % DIR_UKB_CALL+'/ukb_cal_chr{1:22}_v2.bed',
                '--bim=%s' % DIR_UKB_CALL+'/ukb_cal_chr{1:22}_v2.bim',
                '--fam=%s' % fam0,
                '--bgenFile=%s' % DIR_UKB_IMP+'/ukb_imp_chrX_v3.bgen',
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