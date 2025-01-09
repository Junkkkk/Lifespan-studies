import pandas as pd
import numpy as np
from sklearn.model_selection import KFold

pheno_covar = pd.read_csv('pheno_covar.csv')
f=pd.read_csv('completeDataID_393833_all.txt', sep='\t', header=None)

p=pd.merge(f, pheno_covar, left_on=[0], right_on=['eid'])

df = p.sample(frac=1, random_state=42).reset_index(drop=True)

# KFold (Fix random seed)
num_folds = 5
kf = KFold(n_splits=num_folds, shuffle=True, random_state=42)

# Generating fold number
folds = []
for fold_num, (train_index, test_index) in enumerate(kf.split(df)):
    fold_data = df.iloc[train_index]
    fold_data['fold'] = fold_num
    folds.append(fold_data)

# Generating files for each fold
for fold_num, fold_data in enumerate(folds):
    fold_data.drop(columns=['fold'], inplace=True)
    fold_data.to_csv(f'pheno_covar_fold{fold_num+1}.csv', index=False)
    fold_data[[0,1]].to_csv(f'completeDataID_393833_all_fold{fold_num+1}.txt', index=False, sep='\t', header=None)