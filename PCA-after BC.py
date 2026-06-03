# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 11:48:07 2024
@author: jix01
edited and re-run by Amrita Nagasuri 9/25/25
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import seaborn as sns
from statannotations.Annotator import Annotator
from matplotlib.lines import Line2D

# --- output directory ---
OUT_DIR = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/PCA-after BC outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------
# Load data
# ---------------------------
df1 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/pos+neg_proteoBatchCorrPCA output/Neg_ProteoBatchCorrected_data.csv')
df2 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/pos+neg_proteoBatchCorrPCA output/Pos_ProteoBatchCorrected_data.csv')
df3 = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv')

print("Neg (df1) shape:", df1.shape)
print("Pos (df2) shape:", df2.shape)
print("Pheno (df3) shape:", df3.shape)


print("Pheno counts by age:\n", df3["age"].value_counts(dropna=False))
controls_in_pheno = df3.loc[df3["age"]=="C","Sample Name"].str.replace(r"_r[12]$","",regex=True).unique()
cases_in_pheno    = df3.loc[df3["age"]=="P","Sample Name"].str.replace(r"_r[12]$","",regex=True).unique()
print(f"Unique cases in pheno: {len(cases_in_pheno)}")
print(f"Unique controls in pheno: {len(controls_in_pheno)}")


print("✅ Loaded files:")
print("Neg (df1) shape:", df1.shape)
print("Pos (df2) shape:", df2.shape)
print("Pheno (df3) shape:", df3.shape)
print("-"*70)

# Identify sample columns (exclude metadata)
NON_SAMPLE = {'chem_id','Alignment ID','Average Rt(min)','Average Mz','Metabolite name','Adduct type',
              'Post curation result','Fill %','MS/MS assigned','Reference RT','Reference m/z','Formula',
              'Ontology','INCHIKEY','SMILES','Annotation tag (VS1.0)','RT matched','m/z matched',
              'MS/MS matched','Comment','Manually modified for quantification','Manually modified for annotation',
              'Isotope tracking parent ID','Isotope tracking weight number','RT similarity','m/z similarity',
              'Simple dot product','Weighted dot product','Reverse dot product','Matched peaks count',
              'Matched peaks percentage','Total score','S/N average','Spectrum reference file name',
              'MS1 isotopic spectrum','MS/MS spectrum','Category'}

sample_cols1 = [c for c in df1.columns if c not in NON_SAMPLE]
sample_cols2 = [c for c in df2.columns if c not in NON_SAMPLE]
print(f"Sample columns — Neg: {len(sample_cols1)}, Pos: {len(sample_cols2)}")
print("First 10 Neg sample cols:", sample_cols1[:10])
print("First 10 Pos sample cols:", sample_cols2[:10])
print("-"*70)

# --- add _P/_C suffixes from pheno file ---
pheno_map = df3.set_index("Sample Name")["age"].to_dict()

def add_pc_suffix(df, pheno_map):
    renamed = {}
    for c in df.columns:
        if c in pheno_map:
            renamed[c] = f"{c}_{pheno_map[c]}"
    df.rename(columns=renamed, inplace=True)

add_pc_suffix(df1, pheno_map)
add_pc_suffix(df2, pheno_map)

has_c1 = [c for c in df1.columns if c.endswith('_C')]
has_p1 = [c for c in df1.columns if c.endswith('_P')]
has_c2 = [c for c in df2.columns if c.endswith('_C')]
has_p2 = [c for c in df2.columns if c.endswith('_P')]
print(f"After suffix shim — Neg: {len(has_p1)} _P, {len(has_c1)} _C | Pos: {len(has_p2)} _P, {len(has_c2)} _C")
print("-"*70)

# --- concat pos+neg vertically ---
df = pd.concat([df1, df2], axis=0)
print("After concat df shape:", df.shape)

# --- batch/age prefix mapping ---
batch_map = df3.set_index("Sample Name")["batch"].to_dict()
age_map   = df3.set_index("Sample Name")["age"].to_dict()

def base_name(col):
    if col.endswith("_P") or col.endswith("_C"):
        return col[:-2]
    return col

column_mapping = {}
for col in df.columns:
    bname = base_name(col)
    if bname in batch_map and bname in age_map:
        column_mapping[col] = f"{batch_map[bname]}_{age_map[bname]}_{bname}"

sample_cols_all = [c for c in df.columns if c not in NON_SAMPLE]
mapped_count = sum(1 for c in sample_cols_all if c in column_mapping)
print(f"Mapping coverage: {mapped_count}/{len(sample_cols_all)} sample columns will be prefixed with batch/age")
if mapped_count < len(sample_cols_all):
    unmapped = [c for c in sample_cols_all if c not in column_mapping]
    print("⚠️ Unmapped sample columns (showing up to 20):", unmapped[:20])
print("-"*70)

df.rename(columns=column_mapping, inplace=True)
print("After renaming with batch-age prefixes (df):", df.shape)
present_batches = sorted(set([c.split('_',1)[0] for c in df.columns if c[0] in {'1','2','3'}]))
print("Batches detected in column names:", present_batches)
print("-"*70)

# --- set chem_id index ---
df = df.set_index('chem_id')
print("After setting chem_id as index:", df.shape)

# --- batch subsets using strict prefix match ---
loc1 = df.loc[:, df.columns.str.match(r'^1_')]
loc2 = df.loc[:, df.columns.str.match(r'^2_')]
loc3 = df.loc[:, df.columns.str.match(r'^3_')]
print(f"loc1 shape: {loc1.shape} loc2 shape: {loc2.shape} loc3 shape: {loc3.shape}")

# --- improved diagnostics ---
batch_cols = {
    '1': df.columns[df.columns.str.match(r'^1_')].tolist(),
    '2': df.columns[df.columns.str.match(r'^2_')].tolist(),
    '3': df.columns[df.columns.str.match(r'^3_')].tolist(),
}
for b, cols in batch_cols.items():
    print(f"Batch {b}: {len(cols)} columns")
    print(f"  Example cols: {cols[:5]}")

sample_only = [c for c in df.columns if c not in NON_SAMPLE]
missing_batch_cols = [c for c in sample_only if not any(c.startswith(f"{b}_") for b in ['1','2','3'])]
print(f"⚠️ Sample columns not assigned to any batch prefix: {len(missing_batch_cols)}")
if missing_batch_cols:
    print("Example missing sample column names:", missing_batch_cols[:10])
print("-"*70)

# --- concat batches ---
dr = pd.concat([loc1, loc2, loc3], axis=1)
print("After concatenating loc1+2+3 (dr):", dr.shape)

# --- transpose ---
dr = dr.T
print("After transposing dr:", dr.shape)

# --- add group/batch colors ---
dr['color_MC'] = np.where(dr.index.str.contains('_P'), 'red','black')
dr['color_batch'] = np.where(dr.index.str.contains('1_'), '#FFD700',
                     np.where(dr.index.str.contains('2_'),'#2E8B57','#006400'))
print("After assigning color_MC and color_batch (dr):", dr.shape)

# --- diagnostic counts ---
cases_rep   = np.sum(dr.index.str.contains('_P'))
controls_rep= np.sum(dr.index.str.contains('_C'))
base_ids = dr.index.to_series().str.replace('_r1','', regex=False).str.replace('_r2','', regex=False)
uniq_cases    = base_ids[base_ids.str.contains('_P')].nunique()
uniq_controls = base_ids[base_ids.str.contains('_C')].nunique()
print(f"Sample rows in 'dr' (with reps): cases={cases_rep}, controls={controls_rep}, total={len(dr)}")
print(f"Unique biological samples in 'dr': cases={uniq_cases}, controls={uniq_controls}, total={uniq_cases+uniq_controls}")
print("-"*70)

# --- split labels/features ---
dr12 = dr.loc[:, 'color_MC':'color_batch']
dr11 = dr.drop(['color_MC', 'color_batch'], axis=1, errors='ignore')
print("dr11 shape (features):", dr11.shape, "| dr12 shape (labels):", dr12.shape)

# --- NaN summary ---
na_by_sample = dr11.isna().sum(axis=1)
print("NaN summary:\n", na_by_sample.describe())
print("-"*70)

dr11_cleaned = dr11.dropna(axis=1)
print("After dropna on columns (dr11_cleaned):", dr11_cleaned.shape)

# --- Standardize ---
X_std = StandardScaler().fit_transform(dr11_cleaned)
print("After StandardScaler transform:", X_std.shape)

# --- PCA ---
pca = PCA(n_components=20)
principalComponents = pca.fit_transform(X_std)
print("After PCA fit_transform, principalComponents shape:", principalComponents.shape)
print("Explained variance ratio (first 5):", np.round(pca.explained_variance_ratio_[:5], 4))
print("-"*70)

# --- plot explained variance ---
features = range(pca.n_components_)
plt.bar(features, pca.explained_variance_ratio_, color='blue')
plt.xlabel('PCA features')
plt.ylabel('Variance %')
plt.xticks(features)
plt.savefig(OUT_DIR / 'PCA features and the variance explained_after.png', dpi=400, bbox_inches="tight")
plt.show()

# --- save PCA components ---
PCA_components = pd.DataFrame(principalComponents)
print("Created PCA_components:", PCA_components.shape)
PCA_components.columns = ['PC'+ str(col) for col in PCA_components.columns]
PCA_components = PCA_components.reset_index(drop=True)
dr12 = dr12.reset_index(drop=True)
PCA_components = pd.concat([PCA_components, dr12], axis=1)

print("Final PCA_components after concat with labels:", PCA_components.shape)
# === Count plotted points (replicates) per group ===
group_counts = PCA_components['color_MC'].value_counts()
print("\n🟢 Number of PCA dots (replicates) being plotted:")
print(group_counts.rename({'red': 'Loss (red)', 'black': 'Control (black)'}))
print(f"Total PCA points plotted: {len(PCA_components)}")

# (Optional) verify counts by batch as well
batch_counts = PCA_components['color_batch'].value_counts()
print("\n📦 Points per batch color:")
print(batch_counts)


PCA_components.to_csv(OUT_DIR / 'PCA_components.csv', index=False)
print("Saved PCA_components.csv ✅")

print("FINAL DATAFRAME DIMENSIONS (PCA_components):", PCA_components.shape)
# (138, 22)
print("FINAL DATAFRAME DIMENSIONS (PCA_components):", PCA_components.columns) 
# ['PC0', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9',
#        'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16', 'PC17', 'PC18',
#        'PC19', 'color_MC', 'color_batch']

print("Unique color_batch:", PCA_components["color_batch"].unique())
print("Unique color_MC:", PCA_components["color_MC"].unique())
print("-"*70)

# --- scatter by sample type ---
custom_lines = [Line2D([0], [0], color='red', marker='o', lw=0),
                Line2D([0], [0], color='black', marker='o', lw=0)]
fig, ax = plt.subplots()
scatter = plt.scatter(PCA_components['PC0'], PCA_components['PC1'], c=PCA_components['color_MC'])
legend1 = ax.legend(custom_lines, ['Loss', 'Control'],
                    title="Groups", loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()
scatter.figure.savefig(OUT_DIR / 'PC1 and PC2 color-coded by sample type_after.png', dpi=400, bbox_inches="tight")

# --- scatter by batch ---
custom_lines = [Line2D([0], [0], color='#FFD700', marker='o', lw=0),
                Line2D([0], [0], color='#2E8B57', marker='o', lw=0),
                Line2D([0], [0], color='#006400', marker='o', lw=0)]
fig, ax = plt.subplots()
scatter = plt.scatter(PCA_components['PC0'], PCA_components['PC1'], c=PCA_components['color_batch'])
legend1 = ax.legend(custom_lines, ['1', '2','3'],
                    title="Batch", loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()
scatter.figure.savefig(OUT_DIR / 'PC1 and PC2 color-coded by shipment_after.png', dpi=400, bbox_inches="tight")

# --- boxplots ---
df_box = PCA_components.copy()
df_box['color_batch'] = df_box['color_batch'].replace({'#FFD700': '1', '#2E8B57': '2', '#006400': '3'})
x = "color_batch"; y = "PC0"
my_pal = {"1": "#FFD700", "2": "#2E8B57", "3": "#006400"}
present = sorted(df_box['color_batch'].dropna().unique().tolist())
print("Batches present in df_box:", present)
if set(present) != {'1','2','3'}:
    print("⚠️ Warning: Not all three batches present in plotting dataframe.")
plt.figure(figsize=(8, 6))
ax = sns.boxplot(data=df_box, x=x, y=y, order=present, palette={k: my_pal[k] for k in present})
plt.xlabel('Batch', fontsize=14)
plt.ylabel('PC1', fontsize=14)
pairs = []
if '1' in present and '2' in present: pairs.append(("1","2"))
if '2' in present and '3' in present: pairs.append(("2","3"))
if '1' in present and '3' in present: pairs.append(("1","3"))
annotator = Annotator(ax, pairs, data=df_box, x=x, y=y, order=present)
annotator.configure(test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
annotator.apply_and_annotate()
plt.savefig(OUT_DIR / 'boxplot_for_PC1_by_shipment_after.png', dpi=400, bbox_inches="tight")
plt.show()

df_box2 = PCA_components.copy()
df_box2['color_MC'] = df_box2['color_MC'].replace({'red': 'Loss', 'black': 'Control'})
x = "color_MC"; y = "PC0"
my_pal2 = {"Loss": "red", "Control": "black"}
order2 = [c for c in ['Loss','Control'] if c in df_box2['color_MC'].unique()]
print("Groups present in df_box2:", order2)
plt.figure(figsize=(8, 6))
ax = sns.boxplot(data=df_box2, x=x, y=y, order=order2, palette={k: my_pal2[k] for k in order2},
                 medianprops={"color": "white", "linewidth": 1})
pairs2 = [("Loss","Control")] if set(order2) == {'Loss','Control'} else []
annotator = Annotator(ax, pairs2, data=df_box2, x=x, y=y, order=order2)
annotator.configure(test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
annotator.apply_and_annotate()
plt.savefig(OUT_DIR / 'boxplot_for_PC1_by_sample_type_after.png', dpi=400, bbox_inches="tight")
plt.show()
