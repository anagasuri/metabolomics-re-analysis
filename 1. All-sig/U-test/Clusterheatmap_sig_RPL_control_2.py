# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 12:44:11 2025

@author: Ji Xiaowen
"""

# edited and rerun by AMRITA NAGASURI 

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import os  

# out_dir = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/U-test/U-test new/Clusterheatmap_sig_RPL_control_2.py_output"
out_dir = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/final figure 1 heatmaps/font size increased"
os.makedirs(out_dir, exist_ok=True)

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/3. exogenous/output_stats/PosPVFStats2d_all_sig_features.csv')
df2 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/3. exogenous/output_stats/NegPVFStats2d_all_sig_features.csv') 

print(df1.shape)
print(df1.columns)
# (3006, 79)

print(df2.shape)
print(df2.columns)
# (1866, 79)

df = pd.concat([df1, df2], axis=0)

df.columns = df.columns.str.replace(r'0_', '').str.replace(r'1_', '')

df = df.set_index('Alignment ID')

# --------------------------- load phenotype ------------------------------- #
df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')
df1['Sample Name'] = df1['Sample Name'].str.replace(r'_r1', '').str.replace(r'_r2', '')
df1 = df1.drop_duplicates(subset='Sample Name')

df3 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx')
df3['Sample Name'] = df3['Sample Name'].str.replace(r'_r1.d', '').str.replace(r'_r2.d', '')
df3 = df3.drop_duplicates(subset=['Sample Name'])
df3 = df3[df3['Sample Name'].str.contains('BI', na=False)]

# change df1 name as well
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
df1.loc[df1['Sample Name'].str.contains('BI', na=False), 'Sample Name'] = (
    df1['Sample Name'].replace(rename_dict)
)
df1 = df1.drop_duplicates(subset=['Sample Name'])
df1['Sample Name'] = (
    df1['Sample Name']
    .str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
    .str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
    .str.replace(r'CTRL-', 'CTRL', regex=True)
    .str.replace(r'(?<=UC-)(0+)', '', regex=True)
)

# rename df columns name
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
columns_to_rename = [col for col in df.columns if 'BI' in col]
filtered_rename_dict = {col: rename_dict[col] for col in columns_to_rename if col in rename_dict}
df = df.rename(columns=filtered_rename_dict)

# further change UC-CTRL name
df.columns = df.columns.str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
df.columns = df.columns.str.replace(r'UC-0*(\d+)A', r'UC-\1A', regex=True)
df.columns = df.columns.str.replace(r'UC-CTRL-0*(\d+)', r'UC-CTRL\1', regex=True)

print(df.columns.values)

# --------------------------- demographics file ---------------------------- #
df2 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Demographics_EPL_final.xlsx')
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-0', 'UC-', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-CTRL-', 'UC-CTRL', regex=True)
df2 = df2.rename(columns={'StudyID': 'Sample Name'})
df2 = df2.drop_duplicates(subset=['Sample Name'])
df2 = df2[df2['Sample Name'] != 'UC-CTRL10']
df2 = df2[df2['RPL '] != 0]

df4 = pd.merge(df1, df2, on='Sample Name', how='left')

# --------------------------- subset matrix -------------------------------- #
dfH = df.loc[:, df.columns.isin(df2['Sample Name'])]

dfH.columns = (
    dfH.columns
    .str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
    .str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
    .str.replace(r'CTRL-', 'CTRL', regex=True)
    .str.replace(r'(?<=UC-)(0+)', '', regex=True)
)

df4 = df4[df4['Sample Name'].isin(dfH.columns)]

df5 = df4[['Sample Name', 'Case', 'ConceptionMode', 'batch']].set_index('Sample Name')
df5 = df5.loc[dfH.columns, :]
df5 = df5.T
dfH_2 = pd.concat([df5, dfH])

# --------------------------- column colors -------------------------------- #
case_colors = np.where(dfH_2.loc['Case'] == 1, 'orange', 'black')

batch_conditions = [
    dfH_2.loc['batch'] == 1.0,
    dfH_2.loc['batch'] == 2.0,
    dfH_2.loc['batch'] == 3.0
]
batch_choices = ['#FFD700', '#2E8B57', '#006400']
batch_colors = np.select(batch_conditions, batch_choices, default='gray')

color_map = {1: "#984EA3", 2: "#FF7F00", 3: "#377EB8"}
ART_colors = dfH_2.loc['ConceptionMode'].map(color_map)

col_colors = pd.DataFrame({
    'Case': case_colors,
    'Batch': batch_colors,
    'Conception Mode': ART_colors
}, index=dfH.columns)

col_colors.columns = [None] * col_colors.shape[1]

print(dfH.shape)
print(dfH.columns)
# (4872, 56)

dfH.to_csv(os.path.join(out_dir, "dfH_final_matrix_RPL_vs_CONTROL.csv"))

# ============================ SEABORN PART ============================ #

sns.set(font_scale=1.3)

g = sns.clustermap(
    dfH,
    cmap='Blues',
    col_colors=col_colors,
    colors_ratio=(0.025, 0.025),
    standard_scale=0,
    dendrogram_ratio=0.08,
    cbar_pos=(-0.07, .45, .03, .2)
)

# -------- FONT SIZE ADJUSTMENTS (ONLY CHANGE) -------- #

# Heatmap axis labels (samples + features)
g.ax_heatmap.tick_params(axis='x', labelsize=18, rotation=90)
g.ax_heatmap.tick_params(axis='y', labelsize=18)

# X-axis label
g.ax_heatmap.set_xlabel("RPL Samples", fontsize=18, labelpad=12)
g.ax_heatmap.set_ylabel("Significant Chemicals", fontsize=18, labelpad=12)

# Colorbar (blue bar) tick labels
g.ax_cbar.tick_params(labelsize=16)

# Colorbar label (if present)
if g.ax_cbar.get_ylabel():
    g.ax_cbar.set_ylabel(g.ax_cbar.get_ylabel(), fontsize=16)

# Column color strip labels (Case / Batch / Conception Mode)
for label in g.ax_col_colors.get_yticklabels():
    label.set_fontsize(16)

# -------------------------------------------------- #

g.fig.set_size_inches(18, 18)

g.savefig(
    os.path.join(out_dir, 'sig_RPL_control_all_sig_features.png'),
    dpi=300,
    transparent=True
)

# -------------------------- feature order output ------------------------- #
row_order = g.dendrogram_row.reordered_ind
row_labels_ordered = dfH.index[row_order]
pd.Series(row_labels_ordered, name='Feature').to_csv(
    os.path.join(out_dir, 'features_ordered.csv'),
    index=False
)
