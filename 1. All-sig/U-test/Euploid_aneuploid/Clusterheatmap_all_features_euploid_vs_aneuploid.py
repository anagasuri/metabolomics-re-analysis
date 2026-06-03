# -*- coding: utf-8 -*-
"""
ALL FEATURES — Euploid vs Aneuploid (CASES ONLY)

POC_karyotype:
1 = Euploid
0 = Aneuploid
99 = NA (excluded)
"""

# having issues with data alignment 

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os

# ============================
# OUTPUT DIR
# ============================
out_dir = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/final fig S4" 
os.makedirs(out_dir, exist_ok=True)

print("\n=== LOADING FILES ===")

# ============================
# LOAD METABOLOMICS
# ============================
pos = pd.read_csv(
    "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/3. exogenous/neg+pos_removebatcheffect_corrected/Pos_removebatcheffect_corrected_CLEAN.csv"
)
neg = pd.read_csv(
    "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/3. exogenous/neg+pos_removebatcheffect_corrected/Neg_removebatcheffect_corrected_CLEAN.csv"
)

df = pd.concat([pos, neg], axis=0)
df = df.set_index("Alignment ID")

# ============================
# COLLAPSE REPLICATES
# ============================
df.columns = (
    df.columns.astype(str)
    .str.replace("_r1", "", regex=False)
    .str.replace("_r2", "", regex=False)
)

df = df.T.groupby(df.T.index).mean().T

print("Metabolomics shape (features × samples):", df.shape)

# ============================
# LOAD PHENO
# ============================
pheno = pd.read_csv(
    "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv"
)

pheno["Sample Name"] = (
    pheno["Sample Name"]
    .astype(str)
    .str.replace("_r1", "", regex=False)
    .str.replace("_r2", "", regex=False)
)

# remove QC only
pheno = pheno[~pheno["Sample Name"].str.contains("QC|DIWB", na=False)]

# ============================
# BI → UC NAME MAP
# ============================
name_map = pd.read_excel(
    "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx"
)

name_map["Sample Name"] = (
    name_map["Sample Name"]
    .astype(str)
    .str.replace("_r1.d", "", regex=False)
    .str.replace("_r2.d", "", regex=False)
)

name_map = name_map[name_map["Sample Name"].str.contains("BI", na=False)]
rename_dict = dict(zip(name_map["Sample Name"], name_map["Sample Name_2"]))

pheno["Sample Name"] = pheno["Sample Name"].replace(rename_dict)
df = df.rename(columns=rename_dict)

# ============================
# CLEAN UC FORMATTING
# ============================
df.columns = (
    df.columns.astype(str)
    .str.replace("UC-0", "UC-", regex=False)
)

# ============================
# LOAD DEMOGRAPHICS
# ============================
demo = pd.read_excel(
    "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Demographics_EPL_final.xlsx"
)

demo["StudyID"] = (
    demo["StudyID"]
    .astype(str)
    .str.replace("UC-0", "UC-", regex=False)
)

demo = demo.rename(columns={"StudyID": "Sample Name"})

# ============================
# MERGE METADATA (LEFT JOIN)
# ============================
meta = pd.merge(pheno, demo, on="Sample Name", how="left")

# ============================
# BIOLOGICAL FILTER (ONLY THIS)
# ============================
meta = meta[
    (meta["Case"] == 1) &
    (meta["POC_karyotype"].isin([0, 1]))
].copy()

print("\n=== CASE COUNTS ===")
print("Total cases:", meta.shape[0])
print("Euploid (1):", (meta["POC_karyotype"] == 1).sum())
print("Aneuploid (0):", (meta["POC_karyotype"] == 0).sum())

# ============================
# INTERSECT WITH METAB
# ============================
dfH_cols = sorted(set(df.columns).intersection(meta["Sample Name"]))
dfH = df[dfH_cols].copy()

print("\nSamples plotted:", dfH.shape[1])
print(dfH.columns.tolist())

# ============================
# BUILD METADATA MATRIX
# ============================
df_meta = meta.set_index("Sample Name")[["POC_karyotype", "batch"]]
df_meta = df_meta[~df_meta.index.duplicated(keep="first")]
df_meta = df_meta.reindex(dfH.columns)

dfH_annot = pd.concat([df_meta.T, dfH], axis=0)

# ============================
# COLOR BARS
# ============================
karyotype_colors = np.where(
    dfH_annot.loc["POC_karyotype"] == 1,
    "red",      # euploid
    "black"     # aneuploid
)

batch_colors = np.select(
    [
        dfH_annot.loc["batch"] == 1,
        dfH_annot.loc["batch"] == 2,
        dfH_annot.loc["batch"] == 3,
    ],
    ["#FFD700", "#2E8B57", "#006400"],
    default="gray"
)

col_colors = pd.DataFrame(
    {"Karyotype": karyotype_colors, "Batch": batch_colors},
    index=dfH.columns
)
col_colors.columns = [None, None]

# ============================
# HEATMAP
# ============================
sns.set(font_scale=1.0)

g = sns.clustermap(
    dfH,
    cmap="Blues",
    col_colors=col_colors,
    standard_scale=0,
    dendrogram_ratio=0.08,
    colors_ratio=(0.025, 0.05),
    figsize=(20, 20),
    cbar_pos=(-0.10, 0.45, 0.03, 0.25),
)

g.fig.suptitle("All Features — Euploid vs Aneuploid (CASES ONLY)", y=1.02)

out_path = os.path.join(out_dir, "all_features_euploid_aneuploid_CASES_ONLY.png")
g.savefig(out_path, dpi=400, bbox_inches="tight")
plt.show()

print("\nSaved to:", out_path)
