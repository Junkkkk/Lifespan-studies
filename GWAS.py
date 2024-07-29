import os
from multiprocessing import Pool
import pandas as pd

def gwas_chr(i):
    if i == 23:
        i='X'

    os.system(f'plink2 --bfile {result_path}/chr{i} --adjust --linear --pheno {result_path}/pheno.txt --out {result_path}/gwas_chr{i} --ci 0.95')


if __name__ == '__main__':
    num_processes = 23
    with Pool(num_processes) as pool:
        pool.map(gwas_chr, range(1, 24))