# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 12:44:11 2025

@author: Ji Xiaowen
"""
# edited and rerun by AMRITA NAGASURI 
# CONFIRMED DATAFRAMES 

# plots both cases and controls 

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/U-test/RPL_NONRPL/RPL_PosPVFStats2d.csv')
df2 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/U-test/RPL_NONRPL/RPL_NegPVFStats2d.csv')

print(df1.shape)
print(df1.columns)
# (168, 79)

print(df2.shape)
print(df2.columns)
# (118, 79)


df = pd.concat([df1, df2], axis=0)

df.columns = df.columns.str.replace(r'0_', '').str.replace(r'1_', '')

df.columns.values
df = df.set_index('Alignment ID')

df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')
df1['Sample Name'] = df1['Sample Name'].str.replace(r'_r1', '').str.replace(r'_r2', '')
df1 = df1.drop_duplicates(subset= 'Sample Name')


df3 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx')
df3['Sample Name'] = df3['Sample Name'].str.replace(r'_r1.d', '').str.replace(r'_r2.d', '')
df3 = df3.drop_duplicates(subset=['Sample Name'])
df3 = df3[df3['Sample Name'].str.contains('BI', na=False)]


#change df1 name as well
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
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-0', 'UC-', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
df2['StudyID'] = df2['StudyID'].str.replace(r'UC-CTRL-', 'UC-CTRL', regex=True)
df2 = df2.rename(columns = {'StudyID':'Sample Name'})
df2 = df2.drop_duplicates(subset=['Sample Name'])
df2 = df2[df2['RPL '] != 0]


df4 = pd.merge(df1, df2, on = 'Sample Name', how = 'left')

print(df.columns.values)
dfH = df.loc[:, 'UC-20A':'p106640'] # change based on data


dfH.columns = dfH.columns.str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
dfH.columns = dfH.columns.str.replace(r'^.*?(UC-CTRL)', r'\1', regex=True)
dfH.columns = dfH.columns.str.replace(r'CTRL-', 'CTRL', regex=True)
dfH.columns = dfH.columns.str.replace(r'(?<=UC-)(0+)', '', regex=True)



df4 = df4[df4['Sample Name'].isin(dfH.columns)]

df5 = df4[['Sample Name', 'RPL ', 'batch']].set_index('Sample Name')
df5 = df5.loc[dfH.columns, :]
df5 = df5.T
dfH_2 = pd.concat([df5, dfH])

dfH_2 = dfH_2.loc[:, dfH_2.columns.isin(df2['Sample Name'])]



case_colors = np.where(dfH_2.loc['RPL '] == 1, 'orange', 'gray')

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
    'RPL ': case_colors,
    'Batch': batch_colors
})

col_colors.index = dfH_2.columns

col_colors.columns = [None] * col_colors.shape[1]


dfH_2 = dfH_2.drop(index=['RPL ', 'batch'])
sns.set(font_scale=1)
g = sns.clustermap(dfH_2, cmap='Blues', 
                   col_colors=col_colors, 
                   #row_colors=row_colors,
                   colors_ratio=(0.025, 0.025), 
                   standard_scale=0, 
                   dendrogram_ratio=0.08, 
                   cbar_pos=(-0.07, .45, .03, .2))


g.savefig('sig_RPL_control.png', dpi=300, transparent=True)

print("FINAL DATAFRAME DIMS:", dfH.shape) 
# (286, 69), 23 controls and 46 cases 

# extract the order of features
row_order = g.dendrogram_row.reordered_ind
row_labels_ordered = dfH.index[row_order]
row_labels_ordered_series = pd.Series(row_labels_ordered, name='Feature')
row_labels_ordered_series.to_csv('features_ordered.csv', index=False)
