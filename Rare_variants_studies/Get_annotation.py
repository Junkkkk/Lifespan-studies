import pandas as pd
import os
import tqdm
import glob
from multiprocessing import Pool

def process_chr(i):
    wes_file_path='/WES'
    output_path='/UKB'

    os.system(f'plink2 --bfile {wes_file_path}/ukb_wes_chr{i} --geno 0.05 --hwe 1e-6 --make-bed --out {output_path}/WES/ukb_wes_chr{i}_geno_0.05_hwe_1e-6 --threads 2')
    os.system(f'plink2 --bfile {output_path}/WES/ukb_wes_chr{i}_geno_0.05_hwe_1e-6 --max-maf 0.01 --make-bed --out {output_path}/WES/ukb_wes_chr{i}_geno_0.05_hwe_1e-6_maf_u0.01 --threads 2')

def process_data(line):
    columns = line.strip().split()
    title = columns[3]
    chm = columns[0]
    return ' '.join(columns), title, chm

def extract_plink(txt_file):
    wes_path='/UKB/WES'
    out_path='/UKB/WES/split_gene'
    txt_file_path='/UKB/WES/split_gene/txt_files'

    gene=txt_file.split('/')[-1].split('.txt')[0].split('_')[-1]
    i=txt_file.split('/')[-1].split('.txt')[0].split('_')[2].split('chr')[1]

    # if int(i) in [10,11,12,13,14,15,16,17,18,19,20,21,22,1,2,3,4]:
    os.system(f'plink2 --bfile {wes_path}/ukb_wes_chr{i}_geno_0.05_hwe_1e-6_maf_u0.01 --extract range {txt_file_path}/ukb_wes_chr{i}_{gene}.txt --make-pgen --out {out_path}/ukb_wes_chr{i}_{gene}')
    os.system(f'cp {out_path}/ukb_wes_chr{i}_{gene}.pvar {out_path}/ukb_wes_chr{i}_{gene}.vcf')

def annotate_snp_VEP(vcf_file):
    gene=vcf_file.split('/')[-1].split('.vcf')[0].split('_')[-1]
    i=vcf_file.split('/')[-1].split('.vcf')[0].split('_')[2].split('chr')[1]

    os.system(f'/vep -i {vcf_file} \
        -o /UKB/WES/annotated/chr{i}_{gene}.annot \
            --assembly GRCh38 --tab --cache --offline --dir_cache cache\
                --no_stats --force_overwrite --verbose --symbol --pick \
                    --dir_plugins /VEP_plugins/ \
                        --plugin SingleLetterAA \
                            --plugin CADD,/VEP_plugins/whole_genome_SNVs.tsv.gz \
                                --plugin REVEL,file=/VEP_plugins/ref_files/new_tabbed_revel_grch38.tsv.gz \
                                    --plugin LoFtool --plugin LoF,loftee_path:/VEP_plugins/,human_ancestor_fa:/VEP_plugins/ref_files/human_ancestor.fa.gz,phylocsf_data:/VEP_plugins/ref_files/phylocsf_gerp.sql \
                                        --plugin AlphaMissense,file=/VEP_plugins/AlphaMissense_hg38.tsv.gz'
    )

def extract_anntation_info(annot_file):
    save_path='/UKB/WES/annotated/annotated_csv'
    gene=annot_file.split('/')[-1].split('.annot')[0].split('_')[-1]
    i=annot_file.split('/')[-1].split('.annot')[0].split('_')[0]

    os.makedirs(f'/UKB/WES/annotated/annotated_csv/{i}', exist_ok=True)

    df=pd.read_csv(f'{annot_file}', header=None, comment='#', sep='\t')

    cols = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 
    'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 
    'IMPACT', 'DISTANCE', 'STRAND', 'FLAGS', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'HGVSp', 'CADD_PHRED', 'CADD_RAW', 
    'REVEL', 'LoFtool', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info', 'am_class', 'am_pathogenicity']

    df.columns = cols
    df = df.loc[df.SYMBOL ==  gene]
    df = df.loc[((df.Consequence != 'intron_variant') #& (df.Consequence != '3_prime_UTR_variant')
                & (df.Consequence != 'synonymous_variant'))]

    df.index = df.Uploaded_variation
    df_res = pd.DataFrame()
    df_res['SNP'] = list(set(df.index))
    for col in ['Existing_variation', 'Consequence', 'Amino_acids', 'Codons', 'SYMBOL', 'HGVSp', 'IMPACT', 'LoF',  'LoF_filter', 'LoF_flags', 'LoF_info', 'am_class']:
        tab = []
        for snp in df_res.SNP:
            vals = df.loc[snp][col]
            if isinstance(vals, str):
                if col != 'HGVSp':
                    tab.append(vals)
                else:
                    if vals == '-':
                        tab.append('-')
                    else:
                        tab.append(vals.split(':')[1])
            else:
                if col != 'HGVSp':
                    tab.append('|'.join(list(set(vals)-set(['-']))))
                else:
                    hgvsp = list(set(vals)-set(['-']))
                    hgvsp = list(set([h.split(':')[1] for h in hgvsp]))
                    
                    tab.append('|'.join(hgvsp))
        df_res[col] = tab

    for col in ['CADD_PHRED',  'REVEL', 'am_pathogenicity']:
        tab = []
        for snp in df_res.SNP:
            vals = df.loc[snp][col]
            if (isinstance(vals, float) or
                (isinstance(vals, str) and vals != '-')):
                tab.append(float(vals))
            elif 'gnomAD' in col and list(set(vals))[0] != '-':
                if len(list(set(vals))) == 1:
                    tab.append(float(list(set(vals))[0]))
                else:
                    tab.append(float(vals))

            else:
                vals = list(set(vals)-set('-'))
                if len(vals) == 1:
                    tab.append(vals[0])
                elif len(vals) == 0:
                    tab.append(0)
        df_res[col] = tab

    df_res.to_csv(f'{save_path}/{i}/{gene}.csv', index=False, header=True)


if __name__ == '__main__':
    # QC
    num_processes = 22 
    with Pool(num_processes) as pool:
        pool.map(process_chr, 1,23)

    ## Split txt to do plink
    data_file = '/glist-hg38'

    with open(data_file, 'r') as f_in:
        for line in tqdm.tqdm(f_in):
            selected_columns, title, chm = process_data(line)
            output_file = f'/UKB/WES/split_gene/txt_files/ukb_wes_chr{chm}_{title}.txt'

            with open(output_file, 'w') as f_out:
                f_out.write(selected_columns)

    # Split Gene
    txt_files=glob.glob('/UKB/WES/split_gene/txt_files/*.txt')
    
    num_processes = 10
    with Pool(num_processes) as pool:
        pool.map(extract_plink, txt_files)

    # Annotate
    vcf_files=glob.glob('/UKB/WES/split_gene/*.vcf')

    num_processes = 20
    with Pool(num_processes) as pool:
        pool.map(annotate_snp_VEP, vcf_files)

    os.makedirs('/UKB/WES/annotated/annotated_csv', exist_ok=True)

    # Extract Annotation file
    annot_files=glob.glob('/UKB/WES/annotated/*.annot')
    
    num_processes = 20
    with Pool(num_processes) as pool:
        pool.map(extract_anntation_info, annot_files)