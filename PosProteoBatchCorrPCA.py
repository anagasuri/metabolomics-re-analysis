#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 22:31:06 2025

@author: dabrahamsson

edited and re-run by Amrita Nagasuri 9/25/25
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# to output to a specific directory
import os 


df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/pos+neg_imputationV2 output/Pos_Clean_IMPUTED.csv')
df1 = df1.set_index('chem_id')
loc1 = df1.loc[:, 'Alignment ID':'MS/MS spectrum']
loc2 = df1.loc[:, df1.columns.str.contains('BI')]

df2 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')
df2 = df2[df2['Sample Name'].isin(loc2.columns)]
df2 = df2.set_index('Sample Name')

loc2 = loc2.loc[:, loc2.columns.isin(df2.index)]
loc2 = np.log10(loc2)

# Assuming 'Batch' column indicates the batch, and rest are features
batches = df2['batch']
data = loc2.T

# Standardize the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(data)

# Perform PCA to identify batch-related variance
pca = PCA(n_components=10)  # Adjust components as necessary
pca_data = pca.fit_transform(scaled_data)

# Identify components correlated with batch
batch_effects = []
for i in range(pca_data.shape[1]):
    correlation = np.corrcoef(pca_data[:, i], batches)[0, 1]
    if abs(correlation) > 0.1:  # Threshold can be adjusted
        batch_effects.append(i)

# Remove batch effects by subtracting batch-correlated components
pca_data[:, batch_effects] = 0

# Reconstruct the data without batch effects
corrected_data = pca.inverse_transform(pca_data)
corrected_data = scaler.inverse_transform(corrected_data)
corrected_data = 10**(corrected_data)

# Save corrected dataset
corrected_df = pd.DataFrame(corrected_data, columns=data.columns)
corrected_df = corrected_df.T
corrected_df.columns = loc2.columns
corrected_df = pd.concat([loc1, corrected_df], axis=1)

# edited by Amrita to save to a directory 
OUT_DIR = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/pos+neg_proteoBatchCorrPCA output"
os.makedirs(OUT_DIR, exist_ok=True)
out_file = os.path.join(OUT_DIR, "Pos_ProteoBatchCorrected_data.csv")
corrected_df.to_csv(out_file)  # keep same behavior (writes index)
print(f"Saved: {out_file}")
