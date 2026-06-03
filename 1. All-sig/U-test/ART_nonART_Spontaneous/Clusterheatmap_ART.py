# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 12:44:11 2025

@author: Ji Xiaowen
"""
# edited and rerun by AMRITA NAGASURI 
# CONFIRMED DATAFRAMES 

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import os

out_dir = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/U-test/U-test new/Clusterheatmap_ART.py_output"
os.makedirs(out_dir, exist_ok=True)

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/regress output/Pos_removebatcheffect_corrected_NEW.csv')
print(df1.head())
df1 = df1.rename(columns={'name': 'Alignment ID'})
df2 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/regress output/Neg_removebatcheffect_corrected_NEW.csv')
print(df2.head())
df = pd.concat([df1, df2], axis=0)
df.columns.values
print(df.head())

print(df1.shape)
# (12137, 139), makes sense since it contains replicates 
print(df2.shape)
# (8947, 139), makes sense since it contains replicates 

df = df.set_index('Alignment ID')

df.columns = df.columns.str.replace(r'_r[12]$', '', regex=True)
df = df.groupby(df.columns, axis=1).mean()

# match name from DTSC to UCSF + batch info
df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')
df1['Sample Name'] = df1['Sample Name'].str.replace(r'_r1', '').str.replace(r'_r2', '')
df1 = df1.drop_duplicates(subset= 'Sample Name')
df3 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx')
df3['Sample Name'] = df3['Sample Name'].str.replace(r'_r1.d', '').str.replace(r'_r2.d', '')
df3 = df3.drop_duplicates(subset=['Sample Name'])
df3 = df3[df3['Sample Name'].str.contains('BI', na=False)]

# rename df columns name
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
columns_to_rename = [col for col in df.columns if 'BI' in col]
filtered_rename_dict = {col: rename_dict[col] for col in columns_to_rename if col in rename_dict}
df = df.rename(columns=filtered_rename_dict)
# further change UC-Ctrl name 
df.columns = df.columns.str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
df.columns = df.columns.str.replace(r'UC-0*(\d+)A', r'UC-\1A', regex=True)


#change df1 name, df1 is for batch information
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
df1.loc[df1['Sample Name'].str.contains('BI', na=False), 'Sample Name'] = \
    df1['Sample Name'].replace(rename_dict)
df1 = df1.drop_duplicates(subset=['Sample Name'])
df1['Sample Name'] = df1['Sample Name'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'CTRL-', 'CTRL', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'(?<=UC-)(0+)', '', regex=True)

# bioinfo match
df2 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Demographics_EPL_final_NEW.xlsx')
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-0', 'UC-', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-CTRL-', 'UC-CTRL', regex=True)
df2 = df2.rename(columns = {'StudyID':'Sample Name'})
df2 = df2.drop_duplicates(subset=['Sample Name'])
df2 = df2[df2['Sample Name'] != 'UC-CTRL10']

df4 = pd.merge(df1, df2, on = 'Sample Name', how = 'left')
print(df.columns.values)

dfH = df.loc[:, 'UC-20A':'p106640'] # change based on data

dfH.columns = dfH.columns.str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
dfH.columns = dfH.columns.str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
dfH.columns = dfH.columns.str.replace(r'CTRL-', 'CTRL', regex=True)
dfH.columns = dfH.columns.str.replace(r'(?<=UC-)(0+)', '', regex=True)

df4 = df4[df4['Sample Name'].isin(dfH.columns)]

df5 = df4[['Sample Name', 'ConceptionMode', 'batch']].set_index('Sample Name')
df5 = df5.loc[dfH.columns, :]
df5 = df5.T
dfH_2 = pd.concat([df5, dfH])

color_map = {1: "purple", 2: "orange", 3: "teal"}
colors = dfH_2.loc['ConceptionMode'].map(color_map)

# For 'Batch' row color: cool blue for Batch 1, light blue for Batch 2, light orange for Batch 3
batch_conditions = [
    dfH_2.loc['batch'] == 1.0,  
    dfH_2.loc['batch'] == 2.0,
    dfH_2.loc['batch'] == 3.0
]
batch_choices = ['#FFD700', '#2E8B57', '#006400']
batch_colors = np.select(batch_conditions, batch_choices, default='gray')



# Combine all row colors
col_colors = pd.DataFrame({
    'ConceptionMode': colors,
    'Batch': batch_colors
})

col_colors.index = dfH.columns

#col_colors.columns = [None] * col_colors.shape[1]

print("FINAL DATA FRAME DIMENSIONS:", dfH.shape)
# FINAL DATA FRAME DIMENSIONS: (21084, 69)
p_cols = [col for col in dfH.columns if col.startswith('p')]

# Print their names
print("Columns starting with 'p':")
print(p_cols)

# Print how many
print(f"\nNumber of columns starting with 'p': {len(p_cols)}") 

sns.set(font_scale=1)
g = sns.clustermap(dfH, cmap='Blues', 
                   col_colors=col_colors, 
                   #row_colors=row_colors,
                   colors_ratio=(0.025, 0.025), 
                   standard_scale=0, 
                   dendrogram_ratio=0.08, 
                   cbar_pos=(-0.07, .45, .03, .2))


# --- Legend handles ---
conception_legend = [
    Patch(facecolor='purple', label='IVF'),
    Patch(facecolor='orange', label='ICSI'),
    Patch(facecolor='teal', label='Natural')
]

batch_legend = [
    Patch(facecolor='#FFD700', label='Batch 1'),
    Patch(facecolor='#2E8B57', label='Batch 2'),
    Patch(facecolor='#006400', label='Batch 3')
]

# --- Position under the colorbar ---
# We'll stack Conception Mode above Batch
leg1 = g.ax_cbar.legend(
    handles=conception_legend,
    title="Conception Mode",
    loc='upper left',
    bbox_to_anchor=(-0.1, -0.25),  # adjust vertical spacing
    frameon=True
)

leg2 = g.ax_cbar.legend(
    handles=batch_legend,
    title="Batch",
    loc='upper left',
    bbox_to_anchor=(-0.1, -0.55),  # slightly lower
    frameon=True
)

g.ax_cbar.add_artist(leg1)  # ensures both appear

# --- Save + display ---
g.fig.set_size_inches(20, 18)
plt.tight_layout()
g.savefig(
    os.path.join(out_dir, 'Clusterheatmap_ART_legends_left.png'),
    dpi=300,
    bbox_inches='tight',
    pad_inches=0.5,
    transparent=True
)
plt.show()



# === Extract feature order ===
row_order = g.dendrogram_row.reordered_ind
row_labels_ordered = dfH.index[row_order]
row_labels_ordered_series = pd.Series(row_labels_ordered, name='Feature')
row_labels_ordered_series.to_csv(os.path.join(out_dir, 'features_ordered.csv'), index=False)

# g.ax_heatmap.add_artist(leg1)


# g.savefig(os.path.join(out_dir, 'Clusterheatmap_ART.png'), dpi=300, transparent=True)


# # extract the order of features
# row_order = g.dendrogram_row.reordered_ind
# row_labels_ordered = dfH.index[row_order]
# row_labels_ordered_series = pd.Series(row_labels_ordered, name='Feature')
# row_labels_ordered_series.to_csv(os.path.join(out_dir, 'features_ordered.csv'), index=False)

# plt.tight_layout()   # ensures legend and labels don’t get cropped
# plt.show()           # actually displays the heatmap