# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 14:56:12 2025

@author: jix01

edited and re-run by Amrita Nagasuri 9/25/25. 
"""

import pandas as pd
import numpy as np
from scipy import stats

# to output file to directory 
import os 


df = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/pos, neg area/POS_Area_1_2025_01_13_18_46_15_CLEANED.csv')
df['chem_id'] = (
                np.round(df['Average Mz'], 3).astype(str) + '_@_' 
                + np.round(df['Average Rt(min)'], 3).astype(str)
                )

df.columns.values

loc1a = df['chem_id']
loc1b = df.loc[:, 'Alignment ID':'MS/MS spectrum']
loc1 = pd.concat([loc1a, loc1b], axis=1)
loc1 = loc1.set_index('chem_id')

loc2a = df.loc[:, 'BI00700_r1':'Plasma MB-B3_pos']
loc3 = pd.concat([loc1a, loc2a], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 3000] = np.nan

loc3['count'] = loc3.count(axis = 1)
loc3['freq'] = loc3['count']/len(loc3.columns)
loc3['freq%'] = loc3['freq']*100
loc3 = loc3[loc3['freq%'] > 70]
loc3 = loc3.loc[:,  loc3.columns.str.contains('BI')]

loc3L = np.log10(loc3)

def fillNaN_with_unifrand(df):
    lower, upper = 0, df.min()
    a = df.values
    m = np.isnan(a) # mask of NaNs
    mu, sigma = df.min(), df.std()
    a[m] = stats.truncnorm.rvs(
          (lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=m.sum())
    return df


loc3Lmod = loc3L.apply(fillNaN_with_unifrand)
loc3Lmod = 10**loc3Lmod

loc1 = loc1.reset_index()
loc3Lmod = loc3Lmod.reset_index()
df = pd.merge(loc1, loc3Lmod, how='inner')
df = df.set_index('chem_id')
df.to_csv('Pos_Clean.csv')

# --- Save output to the specified directory ---
OUT_DIR = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/pos+neg_imputationV2 output"
os.makedirs(OUT_DIR, exist_ok=True)
out_file = os.path.join(OUT_DIR, "Pos_Clean_IMPUTED.csv")
df.to_csv(out_file)  # unchanged behavior: writes with index
print(f"Saved: {out_file}")