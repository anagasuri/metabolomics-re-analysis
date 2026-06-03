# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 18:25:56 2025

@author: Ji Xiaowen

edited and re-run by Amrita Nagasuri 9/25/25
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# --- added: output directory ---
OUT_DIR = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/PCA+other variables output")
OUT_DIR.mkdir(parents=True, exist_ok=True)
# --------------------------------

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/pos+neg_proteoBatchCorrPCA output/Neg_ProteoBatchCorrected_data.csv')
df2 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/pos+neg_proteoBatchCorrPCA output/Pos_ProteoBatchCorrected_data.csv')
df3 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')

df = pd.concat([df1, df2], axis=0)

column_mapping = {
    col: f"{df3.loc[df3['Sample Name'] == col, 'batch'].values[0]}_{df3.loc[df3['Sample Name'] == col, 'age'].values[0]}_{col}"
    for col in df.columns if col in df3['Sample Name'].values
}

df.rename(columns=column_mapping, inplace=True)

df.columns.values


df = df.set_index('chem_id')
loc1 = df.loc[:, df.columns.str.contains('1_')]
loc2 = df.loc[:, df.columns.str.contains('2_')]
loc3 = df.loc[:, df.columns.str.contains('3_')]
dr = pd.concat([loc1, loc2, loc3], axis=1)
dr.columns = dr.columns.str.replace('_C', '_C')
dr.columns = dr.columns.str.replace('_P', '_P')

dr = dr.T

dr['color_MC'] = np.where(dr.index.str.contains('_P'), 'red','darkgray')
dr['color_batch'] = np.where(dr.index.str.contains('1_'), 'darkgray', 
                    np.where(dr.index.str.contains('2_'),'dodgerblue',
                    'red'))

dr12 = dr.loc[:, 'color_MC':'color_batch']
dr11 = dr.drop(['color_MC', 'color_batch'], axis=1)


# Check for NaNs in dr11
print("NaN counts per column:\n", dr11.isna().sum())

#Remove columns with any NaN values
dr11_cleaned = dr11.dropna(axis=1)

# Standardize the data to have a mean of ~0 and a variance of 1
X_std = StandardScaler().fit_transform(dr11_cleaned)

# Create a PCA instance: pca
pca = PCA(n_components=20)
principalComponents = pca.fit_transform(X_std)

pca.explained_variance_ratio_


# Save components to a DataFrame
PCA_components = pd.DataFrame(principalComponents)



PCA_components.columns = ['PC'+ str(col) for col in PCA_components.columns]

PCA_components.index = dr11_cleaned.index  # Keep the original index
PCA_components = pd.concat([PCA_components, dr12], axis=1)
PCA_components = PCA_components.reset_index()


PCA_components['index'] = PCA_components['index'].str[4:]
PCA_components = PCA_components[PCA_components['index'].str.contains('_r1')]
PCA_components['index'] = PCA_components['index'].str.replace('_r1', '', regex=True)

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')
df1['Sample Name'] = df1['Sample Name'].str.replace(r'_r1', '').str.replace(r'_r2', '')


# bioinfo match
df2 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/metadata/Demographics_EPL_BDYAnalysis.xlsx')
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-0', 'UC-', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-CTRL-', 'UC-CTRL', regex=True)
df2 = df2.rename(columns = {'StudyID':'Sample Name'})
df2 = df2.drop_duplicates(subset=['Sample Name'])


df3 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/metadata/Name Match.xlsx')
df3['Sample Name'] = df3['Sample Name'].str.replace(r'_r1.d', '').str.replace(r'_r2.d', '')
df3 = df3.drop_duplicates(subset=['Sample Name'])
df3 = df3[df3['Sample Name'].str.contains('BI', na=False)]



rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
df1.loc[df1['Sample Name'].str.contains('BI', na=False), 'Sample Name'] = \
    df1['Sample Name'].replace(rename_dict)
df1 = df1.drop_duplicates(subset=['Sample Name'])
df1['Sample Name'] = df1['Sample Name'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'CTRL-', 'CTRL', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'(?<=UC-)(0+)', '', regex=True)


PCA_components['index'] = PCA_components['index'].replace(rename_dict)

PCA_components['index']  = PCA_components['index'] .str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
PCA_components['index']  = PCA_components['index'] .str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
PCA_components['index']  = PCA_components['index'] .str.replace(r'CTRL-', 'CTRL', regex=True)
PCA_components['index']  = PCA_components['index'] .str.replace(r'(?<=UC-)(0+)', '', regex=True)

PCA_components.rename(columns={'index': 'Sample Name'}, inplace=True)

df = pd.merge(PCA_components, df2, on = 'Sample Name', how ='left')
df['SAB'] = df['SAB'].fillna(0)
df = df.dropna(subset=['Age'])

df.to_csv(OUT_DIR / 'PCA_components_+otherVars.csv', index=False)
print(f"Saved: {OUT_DIR / 'PCA_components_+otherVars.csv'}")
