import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os
import sys

# ======================================================
# Setup
# ======================================================
output_dir = '/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/regress output'
os.makedirs(output_dir, exist_ok=True)

# ======================================================
# Load metabolite data
# ======================================================
df = pd.read_csv(os.path.join(output_dir, 'PosPVFStats2c_NEW.csv')).drop(columns=['Unnamed: 0'], errors='ignore')
df.columns = df.columns.str.replace(r'0_', '', regex=True).str.replace(r'1_', '', regex=True)
df = df.set_index('chem_id')
print(df.shape)

# Select BI sample columns
sample_cols = [c for c in df.columns if c.startswith('BI')]
dfH = df[sample_cols]
print(f"\n[Diagnostics] Selected {len(dfH.columns)} raw sample columns (BI IDs).")

# ======================================================
# Load phenotype metadata
# ======================================================
df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')

# ======================================================
# Load name match mapping
# ======================================================
df3 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx')
df3['Sample Name'] = df3['Sample Name'].str.replace(r'(_r1|_r2)?\.d$', '', regex=True)
df3['Sample Name_2'] = df3['Sample Name_2'].str.replace(r'(_r1|_r2)?\.d$', '', regex=True)
rename_dict = dict(zip(df3['Sample Name'], df3['Sample Name_2']))

# Normalize
df1['Sample Name'] = df1['Sample Name'].str.replace(r'(_r1|_r2)?$', '', regex=True)
dfH.columns = dfH.columns.str.replace(r'(_r1|_r2)?$', '', regex=True)

# Apply mapping
df1['Sample Name'] = df1['Sample Name'].replace(rename_dict)
dfH = dfH.rename(columns={k: v for k, v in rename_dict.items() if k in dfH.columns})

print("\n[Diagnostics] After name-matching:")
print(dfH.columns[:20].tolist())
print("Count unique sample names:", len(dfH.columns))

# ======================================================
# Load demographics
# ======================================================
df2 = pd.read_excel('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Demographics_EPL_BDYAnalysis.xlsx')
df2.columns = df2.columns.str.strip()
df2['StudyID'] = df2['StudyID'].astype(str)
df2['StudyID'] = (
    df2['StudyID']
    .str.replace(r'UC-0', 'UC-', regex=True)
    .str.replace(r'(UC-CTRL-)(0+)(\d+)', r'\1\3', regex=True)
    .str.replace(r'UC-CTRL-', 'UC-CTRL', regex=True)
)
df2 = df2.rename(columns={'StudyID': 'Sample Name'})
df2 = df2.drop_duplicates(subset=['Sample Name'])
df2 = df2[df2['Sample Name'] != 'UC-CTRL10']

# ======================================================
# Merge phenotype + demographics
# ======================================================
df4 = pd.merge(df1, df2, on='Sample Name', how='left')
df4.columns = df4.columns.str.strip()

# Merge RPL fields if duplicated
if {'RPL_x', 'RPL_y'}.issubset(df4.columns):
    df4['RPL'] = df4['RPL_y'].combine_first(df4['RPL_x'])
    df4.drop(columns=['RPL_x', 'RPL_y'], inplace=True)
elif 'RPL_y' in df4.columns:
    df4.rename(columns={'RPL_y': 'RPL'}, inplace=True)
elif 'RPL_x' in df4.columns:
    df4.rename(columns={'RPL_x': 'RPL'}, inplace=True)

print("\n[Diagnostics] Columns in merged metadata:", df4.columns.tolist())

# ======================================================
# Overlap check
# ======================================================
overlap = set(dfH.columns).intersection(set(df4['Sample Name']))
print(f"\n[Diagnostics] Overlapping samples: {len(overlap)}")

if len(overlap) == 0:
    print("⚠️ No overlapping names found!")
    print("Example dfH columns:", list(dfH.columns[:10]))
    print("Example df4 Sample Names:", list(df4['Sample Name'].unique()[:10]))
    sys.exit("[Exit] Please check naming consistency.")

dfH = dfH[list(overlap)]
df4 = df4[df4['Sample Name'].isin(overlap)]

# ======================================================
# Handle RPL
# ======================================================
df4['RPL'] = df4.apply(
    lambda row: 0 if pd.isna(row['RPL']) and 'CTRL' in str(row['Sample Name']).upper()
    else (2 if pd.isna(row['RPL']) else row['RPL']),
    axis=1
)

# ======================================================
# Prepare metadata matrix
# ======================================================
df5 = df4[['Sample Name', 'Case', 'RPL', 'batch']].set_index('Sample Name')

# Ensure unique index names
df5 = df5[~df5.index.duplicated(keep='first')]
dfH = dfH.loc[:, ~dfH.columns.duplicated()]

# Align both before concat
df5 = df5.loc[dfH.columns, :]
df5 = df5.T

# Ensure df5 index names are unique
df5.index = pd.Index(df5.index.map(str))
dfH.index = pd.Index(dfH.index.map(str))

# Concatenate metadata + data
dfH_2 = pd.concat([df5, dfH], axis=0)

# ======================================================
# Fill missing Case values
# ======================================================
for col in dfH_2.columns:
    if pd.isna(dfH_2.loc['Case', col]):
        if 'CTRL' in col.upper():
            dfH_2.loc['Case', col] = 0
        else:
            dfH_2.loc['Case', col] = 1

# ======================================================
# Color bars
# ======================================================
case_colors = np.where(dfH_2.loc['Case'] == 1, 'red', 'black')

batch_conditions = [
    dfH_2.loc['batch'] == 1.0,
    dfH_2.loc['batch'] == 2.0,
    dfH_2.loc['batch'] == 3.0
]
batch_choices = ['#FFD700', '#2E8B57', '#006400']
batch_colors = np.select(batch_conditions, batch_choices, default='gray')

rpl_conditions = [
    dfH_2.loc['RPL'] == 0,
    dfH_2.loc['RPL'] == 1,
]
rpl_choices = ['pink', 'purple']
rpl_colors = np.select(rpl_conditions, rpl_choices, default='gray')

col_colors = pd.DataFrame({
    'Case': case_colors,
    'Batch': batch_colors,
    'RPL': rpl_colors
})
col_colors.index = dfH.columns
col_colors.columns = [None] * col_colors.shape[1]

# ======================================================
# Heatmap
# ======================================================
sns.set(font_scale=1)

print("\n[Diagnostics] Final dfH shape used for plotting:", dfH.shape)
# (6194, 69)
print(f"→ Features (rows): {dfH.shape[0]} | Samples (columns): {dfH.shape[1]}")
print("Unique sample names in dfH:", len(dfH.columns.unique()))
print("First few sample names:", dfH.columns[:10].tolist())

g = sns.clustermap(
    dfH,
    cmap='Blues',
    col_colors=col_colors,
    figsize=(12, 50),  # wider & taller
    colors_ratio=(0.025, 0.025),
    standard_scale=0,
    dendrogram_ratio=0.08,
    cbar_pos=(-0.07, .45, .03, .2)
)

print("Heatmap input non-NaN rows:", dfH.dropna(how='all').shape[0])


output_path = os.path.join(output_dir, 'sig_case_control_2.png')
g.savefig(output_path, dpi=300, transparent=True)
print(f"\n[Done] Figure saved → {output_path}")
