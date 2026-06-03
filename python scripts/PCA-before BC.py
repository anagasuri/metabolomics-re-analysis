# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 11:48:07 2024

@author: jix01
edited by: Amrita
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path

# --- Output directory ---
OUT = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC")
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------
# Load input files
# ---------------------------
df1 = pd.read_csv('/Users/amritanagasuri/Downloads/Q-TOF Data Analysis/2. imput-PCA/Before BC/neg,pos,pheno_inputs/Neg_Clean_NEWWITHREPS.csv')
df2 = pd.read_csv('/Users/amritanagasuri/Downloads/Q-TOF Data Analysis/2. imput-PCA/Before BC/neg,pos,pheno_inputs/Pos_Clean_NEWWITHREPS.csv')
df3 = pd.read_csv('/Users/amritanagasuri/Downloads/Q-TOF Data Analysis/2. imput-PCA/Before BC/neg,pos,pheno_inputs/pheno-clean-NEWWITHREPS.csv')

print("✅ Loaded files:")
print("Neg (df1):", df1.shape)
print("Pos (df2):", df2.shape)
print("Pheno (df3):", df3.shape)
print("-" * 70)

# ---------------------------
# Concatenate positive and negative data
# ---------------------------
df = pd.concat([df1, df2], axis=0)
print("After concat df shape:", df.shape)
print("-" * 70)

# ---------------------------
# Build mapping for batch + phenotype
# ---------------------------
column_mapping = {}
for col in df.columns:
    if col in df3['Sample Name'].values:
        match = df3.loc[df3['Sample Name'] == col]
        batch = match['batch'].values[0]
        age = match['age'].values[0]
        column_mapping[col] = f"{batch}_{age}_{col}"

mapped_count = len(column_mapping)
print(f"Mapped {mapped_count} / {len(df.columns)} columns to pheno-clean entries.")
if mapped_count < len(df.columns):
    missing = [c for c in df.columns if c not in column_mapping]
    print("⚠️ Unmapped columns (first 10):", missing[:10])
print("-" * 70)

df.rename(columns=column_mapping, inplace=True)
print("After renaming with batch-age prefixes:", df.shape)

# ---------------------------
# Separate by batch number
# ---------------------------
if 'chem_id' in df.columns:
    df = df.set_index('chem_id')
print("After setting chem_id as index:", df.shape)

loc1 = df.loc[:, df.columns.str.contains('1_')]
loc2 = df.loc[:, df.columns.str.contains('2_')]
loc3 = df.loc[:, df.columns.str.contains('3_')]

print(f"loc1 shape: {loc1.shape} | loc2 shape: {loc2.shape} | loc3 shape: {loc3.shape}")
print("Batch 1 example cols:", loc1.columns[:5].tolist())
print("Batch 2 example cols:", loc2.columns[:5].tolist())
print("Batch 3 example cols:", loc3.columns[:5].tolist())
print("-" * 70)

# Combine all batches
dr = pd.concat([loc1, loc2, loc3], axis=1)
print("After concatenating loc1+2+3 (dr):", dr.shape)

# Replace suffixes (no-op but retained)
dr.columns = dr.columns.str.replace('_C', '_C', regex=False)
dr.columns = dr.columns.str.replace('_P', '_P', regex=False)

# Transpose to have samples as rows
dr = dr.T
print("After transpose (dr):", dr.shape)

# ---------------------------
# Add color codes
# ---------------------------
dr['color_MC'] = np.where(dr.index.str.contains('_P'), 'red', 'black')
dr['color_batch'] = np.where(dr.index.str.contains('1_'), '#FFD700',
                     np.where(dr.index.str.contains('2_'), '#2E8B57',
                     '#006400'))
print("After assigning color_MC & color_batch:", dr.shape)

# Sample counts
cases_rep = np.sum(dr.index.str.contains('_P'))
controls_rep = np.sum(dr.index.str.contains('_C'))
base_ids = dr.index.to_series().str.replace('_r1','', regex=False).str.replace('_r2','', regex=False)
uniq_cases = base_ids[base_ids.str.contains('_P')].nunique()
uniq_controls = base_ids[base_ids.str.contains('_C')].nunique()

print(f"Samples (with reps): cases={cases_rep}, controls={controls_rep}, total={len(dr)}")
print(f"Unique biological samples: cases={uniq_cases}, controls={uniq_controls}, total={uniq_cases+uniq_controls}")
print("-" * 70)

# ---------------------------
# Split into data and labels
# ---------------------------
dr12 = dr.loc[:, 'color_MC':'color_batch']
dr11 = dr.drop(['color_MC', 'color_batch'], axis=1)

print("dr11 (features):", dr11.shape)
print("dr12 (labels):", dr12.shape)
print("-" * 70)

# ---------------------------
# Handle missing values
# ---------------------------
print("NaN summary before dropna:")
print(dr11.isna().sum().describe())
dr11_cleaned = dr11.dropna(axis=1)
print("After dropna (dr11_cleaned):", dr11_cleaned.shape)
print("-" * 70)

# ---------------------------
# Standardize & PCA
# ---------------------------
X_std = StandardScaler().fit_transform(dr11_cleaned)
print("After StandardScaler:", X_std.shape)

pca = PCA(n_components=20)
principalComponents = pca.fit_transform(X_std)
print("After PCA fit_transform:", principalComponents.shape)
print("Explained variance ratio (first 5):", np.round(pca.explained_variance_ratio_[:5], 4))
print("-" * 70)

# ---------------------------
# Create PCA components DataFrame
# ---------------------------
PCA_components = pd.DataFrame(principalComponents)
print("Created PCA_components:", PCA_components.shape)

# Add color metadata
PCA_components.columns = ['PC' + str(i) for i in PCA_components.columns]
PCA_components = PCA_components.reset_index(drop=True)
dr12 = dr12.reset_index(drop=True)
PCA_components = pd.concat([PCA_components, dr12], axis=1)
print("Final PCA_components with color columns:", PCA_components.shape)
PCA_components.to_csv(OUT / "PCA_components_beforeBC.csv", index=False)
print("Saved PCA_components_beforeBC.csv ✅")
print("-" * 70)

# ---------------------------
# Visual diagnostics (kept original)
# ---------------------------
plt.bar(range(pca.n_components_), pca.explained_variance_ratio_, color='blue')
plt.xlabel('PCA features')
plt.ylabel('Variance %')
plt.savefig(OUT / "PCA_features_variance_explained_beforeBC.png", dpi=400, bbox_inches="tight")
plt.close()

plt.scatter(PCA_components['PC0'], PCA_components['PC1'], alpha=.4, color='blue')
plt.xlabel('PC1'); plt.ylabel('PC2')
plt.savefig(OUT / "PCA_scatter_beforeBC.png", dpi=400, bbox_inches="tight")
plt.close()

print("✅ PCA-before BC completed successfully.")
