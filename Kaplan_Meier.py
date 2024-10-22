import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from matplotlib import pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter, statistics

pheno_covar=pd.read_csv('pheno_covar.csv')

res=pd.read_csv(f'_all_results.csv')
res=res[res.P<7.4e-07]

genes_variant={'gene': res['gene'], 'variant_type': res['SNPs_category']}

for gene, variant_type in tqdm(zip(genes_variant['gene'], genes_variant['variant_type'])):
    chm='chr'+str(res.loc[res['gene']==gene,'chr'].iloc[0])
    
    raw_file=pd.read_csv(f'{chm}_{gene}/ukb_wes_{chm}_{gene}_{variant_type}.raw', sep=' ')
    row_sums = raw_file.iloc[:, 6:].sum(axis=1)
    raw_file[f'{variant_type}_carrier'] = (row_sums != 0).astype(int)

    df=pd.merge(pheno_covar, raw_file, on=['FID','IID'])
    
    data_0 = df[df[f"{variant_type}_carrier"] == 0]
    data_1 = df[df[f"{variant_type}_carrier"] == 1]
    
    kmf_0 = KaplanMeierFitter()
    kmf_1 = KaplanMeierFitter()
    
    kmf_0.fit(data_0['last_known_age'], event_observed=data_0['Death'], entry=data_0['Age_1st_visit'], label=f"Non Carrier")
    kmf_1.fit(data_1['last_known_age'], event_observed=data_1['Death'], entry=data_1['Age_1st_visit'], label=f"Carrier")
    
    plt.figure(figsize=(8, 8))
    
    color_blue = '#7ECBA5'
    color_red = '#d98196'
    
    kmf_0.plot(ci_show=True, color=color_blue, loc=slice(40, None))
    kmf_1.plot(ci_show=True, color=color_red, loc=slice(40, None))
    
    plt.xticks([40, 50, 60, 70, 80, 90])
    plt.yticks([0, 0.25, 0.5, 0.75, 1])

    plt.legend(loc='upper left', bbox_to_anchor=(0, 0.12), frameon=False, ncol=2, prop={'size': 14})
    plt.title(f'{gene}', loc='center', pad=10, fontweight='bold', fontsize=14)
    plt.xlabel('Age', fontsize=14)
    plt.ylabel('Survival Probability', fontsize=12)

    plt.tick_params(axis='x', labelsize=10)
    plt.tick_params(axis='y', labelsize=10)

    plt.savefig(f'{save_path}/{variant_type}_{gene}.pdf', format='pdf')
    plt.show()
    plt.close()