
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 12:44:11 2025

@author: Ji Xiaowen
"""
# CONFIRMED DATAFRAMES 

import numpy as np
import pandas as pd 
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import os 

# edited and rerun by AMRITA NAGASURI 
# CONFIRMED DATAFRAMES 

# only plot cases in this heatmap 

out_dir = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/U-test/U-test new/Clusterheatmap_sig_euploid_aneuploid.py_output"
os.makedirs(out_dir, exist_ok=True)

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/U-test/U-test new/Stats-U-Test.py_output/Aneuploid_PosPVFStats2d.csv')
df2 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/U-test/U-test new/Stats-U-Test.py_output/Aneuploid_NegPVFStats2d.csv')

print(df1.shape)
# (57, 79), 10 metadata cols 
print(df2.shape)
# (68, 79), 10 metadata cols 

df = pd.concat([df1, df2], axis=0)
print(df.shape)
# (125, 80), 1 extra column because a new column 'Alignment ID.1' was created during concatenations 

df.columns = df.columns.str.replace(r'0_', '').str.replace(r'1_', '')

df.columns.values
# first col is named name and not alignment id 
df = df.set_index('name')

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')
df1['Sample Name'] = df1['Sample Name'].str.replace(r'_r1', '').str.replace(r'_r2', '')
df1 = df1.drop_duplicates(subset= 'Sample Name')


df3 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx')
df3['Sample Name'] = df3['Sample Name'].str.replace(r'_r1.d', '').str.replace(r'_r2.d', '')
df3 = df3.drop_duplicates(subset=['Sample Name'])
df3 = df3[df3['Sample Name'].str.contains('BI', na=False)]

#change df1 name as well based on df3
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
df1.loc[df1['Sample Name'].str.contains('BI', na=False), 'Sample Name'] = \
    df1['Sample Name'].replace(rename_dict)
df1 = df1.drop_duplicates(subset=['Sample Name'])
df1['Sample Name'] = df1['Sample Name'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'CTRL-', 'CTRL', regex=True)
df1['Sample Name'] = df1['Sample Name'].str.replace(r'(?<=UC-)(0+)', '', regex=True)


# rename df columns name
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))
columns_to_rename = [col for col in df.columns if 'BI' in col]
filtered_rename_dict = {col: rename_dict[col] for col in columns_to_rename if col in rename_dict}
df = df.rename(columns=filtered_rename_dict)

# further change UC-Ctrl name 
df.columns = df.columns.str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)



# bioinfo match
df2 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Demographics_EPL_final.xlsx')
print(df2.columns)
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-0', 'UC-', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-CTRL-', 'UC-CTRL', regex=True)
df2 = df2.rename(columns = {'StudyID':'Sample Name'})
df2 = df2.drop_duplicates(subset=['Sample Name'])
df2 = df2[df2['Case'] == 1]
df2 = df2[df2['POC_karyotype'] != 99]

df4 = pd.merge(df2, df1, on = 'Sample Name', how = 'left')

print(df.columns.values)
dfH = df.loc[:, 'UC-20A':'p106640'] # change based on data

dfH.columns = dfH.columns.str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
dfH.columns = dfH.columns.str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
dfH.columns = dfH.columns.str.replace(r'CTRL-', 'CTRL', regex=True)
dfH.columns = dfH.columns.str.replace(r'(?<=UC-)(0+)', '', regex=True)

# --- Normalize sample names to consistent zero-padding ---
# --- Normalize sample names to consistent zero-padding ---
# --- Normalize sample names to consistent zero-padding ---

# Fix dfH columns (metabolomics samples)
dfH.columns = dfH.columns.str.replace(r'UC-(\d{2})(?=A)', r'UC-0\1', regex=True)
dfH.columns = dfH.columns.str.replace(r'UC-CTRL-(\d{2,3})', 
                                      lambda m: f"UC-CTRL-{int(m.group(1)):05d}", regex=True)

# Fix df4 sample names (metadata)
df4['Sample Name'] = df4['Sample Name'].str.replace(r'UC-(\d{2})(?=A)', r'UC-0\1', regex=True)
df4['Sample Name'] = df4['Sample Name'].str.replace(r'UC-CTRL-(\d{2,3})', 
                                                   lambda m: f"UC-CTRL-{int(m.group(1)):05d}", regex=True)

# --- Normalize sample names to consistent zero-padding ---
# --- Normalize sample names to consistent zero-padding ---
# --- Normalize sample names to consistent zero-padding ---

dfH = dfH.loc[:, dfH.columns.isin(df4['Sample Name'])]

df5 = df4[['Sample Name', 'POC_karyotype', 'batch']].set_index('Sample Name')
df5 = df5.loc[dfH.columns, :]
df5 = df5.T
dfH_2 = pd.concat([df5, dfH])


case_colors = np.where(dfH_2.loc['POC_karyotype'] == 0, 'blue', 'mediumpurple')

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
    'RPL': case_colors,
    'Batch': batch_colors
})

col_colors.index = dfH.columns

col_colors.columns = [None] * col_colors.shape[1]

print("FINAL DATAFRAME DIMS:", dfH.shape)
# (125, 39), there were 7 unknown, untested case samples 
print("FINAL DATAFRAME COLS:", dfH.columns)
# ['UC-020A', 'UC-023A', 'UC-024A', 'UC-028A', 'UC-034A', 'UC-035A',
#        'UC-038A', 'UC-041A', 'UC-044A', 'UC-046A', 'UC-048A', 'UC-049A',
#        'UC-050A', 'UC-051A', 'UC-053A', 'UC-054A', 'UC-055A', 'UC-056A',
#        'UC-061A', 'UC-062A', 'UC-064A', 'UC-066A', 'UC-070A', 'UC-071A',
#        'UC-072A', 'UC-075A', 'UC-076A', 'UC-080A', 'UC-082A', 'UC-085A',
#        'UC-086A', 'UC-087A', 'UC-092A', 'UC-097A', 'UC-098A', 'UC-099A',
#        'UC-101A', 'UC-103A', 'UC-106A']


sns.set(font_scale=1)
g = sns.clustermap(dfH, cmap='Blues', 
                   col_colors=col_colors, 
                   #row_colors=row_colors,
                   colors_ratio=(0.025, 0.025), 
                   standard_scale=0, 
                   dendrogram_ratio=0.08, 
                   cbar_pos=(-0.07, .45, .03, .2))


# --- Legend handles ---
karyotype_legend = [
    Patch(facecolor='blue', label='Euploid'),
    Patch(facecolor='mediumpurple', label='Aneuploid')
]

batch_legend = [
    Patch(facecolor='#FFD700', label='Batch 1'),
    Patch(facecolor='#2E8B57', label='Batch 2'),
    Patch(facecolor='#006400', label='Batch 3')
]

# --- Position under the colorbar ---
leg1 = g.ax_cbar.legend(
    handles=karyotype_legend,
    title="Karyotype (Top Bar)",
    loc='upper left',
    bbox_to_anchor=(-0.1, -0.25),  # adjust spacing as needed
    frameon=True
)

leg2 = g.ax_cbar.legend(
    handles=batch_legend,
    title="Batch (Second Bar)",
    loc='upper left',
    bbox_to_anchor=(-0.1, -0.55),
    frameon=True
)

g.ax_cbar.add_artist(leg1)

# --- Save + display ---
g.fig.set_size_inches(20, 18)
plt.tight_layout()
g.savefig(
    os.path.join(out_dir, 'sig_euploid_aneuploid_legends_left.png'),
    dpi=300,
    bbox_inches='tight',
    pad_inches=0.5,
    transparent=True
)
plt.show()


g.savefig(os.path.join(out_dir, 'sig_euploid_aneuploid.png'), dpi=300, transparent=True)




