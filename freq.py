# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:57:20 2024

@author: Abrahd05
"""
# re-run/modifed input files by AMRITA NAGASURI for re-analysis

import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

# === Add output directory ===
OUT_DIR = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/3. Frequency/freq.py output"
os.makedirs(OUT_DIR, exist_ok=True)

df = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/3. Frequency/CLEANED_final_annotated_exogenous.csv')
print(df.shape)
# (687, 79) and 10 cols are metadata, so 69 sample cols = 46 cases, 23 controls, file does not contain duplicates 
# df = df.set_index('Compound_id') --> changed to Alignment ID
df = df.set_index('Alignment ID')

print("dimensions after reading in:", df.shape)
# dimensions after reading in: (687, 78) 
# because alignmnet id is now index,removing it as a col, so 9 metadata cols 

# make the average between r1 and r2, and substract MB blanks
df1 =  df.loc[:, df.columns.str.contains('BI')]
df1[df1 < 500] = np.nan

print("dimensions of df1:", df1.shape)
# dimensions of df1: (687, 69)

#rename the columns for control -pos
df_pheno = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')
print(df_pheno.shape)
df_pheno['Sample Name'] = df_pheno['Sample Name'].str.replace('_r1', '').str.replace('_r2', '')
print(df_pheno.shape)
df_pheno = df_pheno.drop_duplicates(subset=['Sample Name'])
print(df_pheno.shape)
column_mapping = {
    col: f"{df_pheno.loc[df_pheno['Sample Name'] == col, 'age'].values[0]}_{col}"
    for col in df1.columns if col in df_pheno['Sample Name'].values
}

df1.rename(columns=column_mapping, inplace=True)

print("dimensions of df1:", df1.shape)

#case
df1_loss =df1.loc[:, df1.columns.str.contains('P_')]
df1_loss['count_loss'] = df1_loss.count(axis = 1)
# adds 1 column
df1_loss['freq_loss'] = df1_loss['count_loss']/len(df1_loss.columns)
# adds 1 column 
df1_loss['freq%_loss'] = df1_loss['freq_loss']*100
# adds 1 column 

print("df1_loss dims:", df1_loss.shape)
# 46 + 3 new cols is 49
print(df1.columns[df1.columns.str.contains('P_')])
# this only looks at cols starting with P_ which is 46 
print(len(df1.columns[df1.columns.str.contains('P_')]))
# prints 46 which is correct! 



#control
df1_control = df1.loc[:, df1.columns.str.contains('C_')]
df1_control['count_control'] = df1_control.count(axis = 1)
# adds 1 column
df1_control['freq_control'] = df1_control['count_control']/len(df1_control.columns)
# adds 1 column
df1_control['freq%_control'] = df1_control['freq_control']*100
# adds 1 column

print("df1_control dims:", df1_control.shape)
# shows (687, 26) with 3 extra cols but its due to above
# real controls count is 23 

# 46 cases + 23 controls = 69 samples total

df2 = pd.concat([df1_loss, df1_control], axis=1)
df2a = df2[['freq%_loss', 'freq%_control']]
df2a.columns = ['Loss (%)', 'Control (%)']


import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(2, 30))
ax = sns.heatmap(
    df2a,
    annot=True,
    fmt=".0f",
    cmap="Blues",
    linewidths=0.5,
    cbar_kws={"shrink": 0.3}
)

colorbar = ax.collections[0].colorbar
colorbar.set_label("Detection Frequency (%)")

plt.xlabel("")  # hide x-axis label if needed
plt.tight_layout()

# === Save output to your directory ===
plt.savefig(os.path.join(OUT_DIR, "freq_heatmap.png"), dpi=300, bbox_inches='tight', transparent=True)
df2.to_csv(os.path.join(OUT_DIR, "detection_frequency_full.csv"))
df2a.to_csv(os.path.join(OUT_DIR, "detection_frequency_summary.csv"))

plt.show()



