# -*- coding: utf-8 -*-
"""
Before Batch Correction — Correlation tables ONLY (NO PLOTS)

Adds StudyID + SAB + RPL by aligning PCA rows with pheno-clean + Name Match + demographics.

Outputs:
- correlation_r_beforeBC_with_SAB_RPL.csv
- correlation_p_beforeBC_with_SAB_RPL.csv
- input_used_beforeBC_with_SAB_RPL.csv
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from pathlib import Path

# -----------------------------
# Paths
# -----------------------------
PCA_FILE   = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/PCA_components_beforeBC.csv")
PHENO_FILE = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv")
NAME_FILE  = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx")
DEMO_FILE  = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Demographics_EPL_final.xlsx")

OUT_DIR = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/correlation output")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Load files
# -----------------------------
pca = pd.read_csv(PCA_FILE)
pheno = pd.read_csv(PHENO_FILE)
name_map = pd.read_excel(NAME_FILE)
demo = pd.read_excel(DEMO_FILE)

# -----------------------------
# Clean column names
# -----------------------------
for df in [pca, pheno, name_map, demo]:
    df.columns = df.columns.str.replace("\u00A0", " ").str.strip()

# Fix demographics typo
if "RPL " in demo.columns and "RPL" not in demo.columns:
    demo = demo.rename(columns={"RPL ": "RPL"})

# Drop junk index column
pca = pca.drop(columns=["Unnamed: 0"], errors="ignore")

# -----------------------------
# Attach Sample Name to PCA by ROW ORDER
# -----------------------------
if "Sample Name" not in pheno.columns:
    raise ValueError(f"pheno file missing 'Sample Name'. Found: {pheno.columns.tolist()}")

n = min(len(pca), len(pheno))
pca = pca.reset_index(drop=True)
pheno = pheno.reset_index(drop=True)

pca.loc[:n-1, "Sample Name"] = pheno.loc[:n-1, "Sample Name"].astype(str)

# -----------------------------
# Map Sample Name -> StudyID
# -----------------------------
required_nm = {"Sample Name", "Sample Name_2"}
if not required_nm.issubset(name_map.columns):
    raise ValueError(f"Name Match missing columns {required_nm}")

pca["SampleName_d"] = pca["Sample Name"].astype(str) + ".d"

name_dict = dict(
    zip(
        name_map["Sample Name"].astype(str),
        name_map["Sample Name_2"].astype(str)
    )
)

pca["StudyID"] = pca["SampleName_d"].map(name_dict).fillna(pca["SampleName_d"])

# -----------------------------
# Merge demographics
# -----------------------------
need_demo = ["StudyID", "SAB", "RPL"]
missing = [c for c in need_demo if c not in demo.columns]
if missing:
    raise ValueError(f"Demographics missing columns: {missing}")

pca = pca.merge(demo[need_demo], on="StudyID", how="left")

# -----------------------------
# Define Sample Type (LOSS vs CONTROL)
# -----------------------------
# 1 = Loss (RPL), 0 = Control
pca["Sample Type"] = pca["RPL"]

print("\nSample Type counts BEFORE filtering:")
print(pca["Sample Type"].value_counts(dropna=False))

# -----------------------------
# Encode Batch ONLY (optional)
# -----------------------------
def factorize_series(s: pd.Series) -> pd.Series:
    return pd.factorize(s.astype(str))[0] + 1

if "color_batch" in pca.columns:
    pca["Batch"] = factorize_series(pca["color_batch"])

# Ensure numeric
for col in ["Sample Type", "Batch", "RPL", "SAB"]:
    if col in pca.columns:
        pca[col] = pd.to_numeric(pca[col], errors="coerce")

# -----------------------------
# Select columns for correlation
# -----------------------------
pc_cols = [c for c in pca.columns if c.startswith("PC")]

df_use = pca[
    pc_cols + ["Sample Type", "Batch", "RPL", "SAB"]
].copy()

# -----------------------------
# Rename PCs: PC0 -> PC1, PC1 -> PC2, ...
# -----------------------------
df_use = df_use.rename(
    columns={f"PC{i}": f"PC{i+1}" for i in range(200)}
)

# Recompute PC columns AFTER renaming
pc_cols = [c for c in df_use.columns if c.startswith("PC")]

# -----------------------------
# Drop rows missing PCs ONLY
# -----------------------------
df_use = df_use.dropna(subset=pc_cols)

# -----------------------------
# Drop samples without Sample Type (recommended)
# -----------------------------
df_use = df_use.dropna(subset=["Sample Type"])

print("\nSample Type counts USED in correlation:")
print(df_use["Sample Type"].value_counts())

assert df_use["Sample Type"].nunique() > 1, \
    "ERROR: Sample Type has no variation after filtering!"

print("\n=== Rows used for correlation ===")
print("N rows used:", len(df_use))

# -----------------------------
# Compute correlation matrix + p-values
# -----------------------------
cols = df_use.columns.tolist()

corr = pd.DataFrame(index=cols, columns=cols, dtype=float)
pvals = pd.DataFrame(index=cols, columns=cols, dtype=float)

for c1 in cols:
    for c2 in cols:
        r, p = pearsonr(df_use[c1], df_use[c2])
        corr.loc[c1, c2] = r
        pvals.loc[c1, c2] = p

# -----------------------------
# Save outputs (NO PLOTS)
# -----------------------------
corr.to_csv(OUT_DIR / "correlation_r_beforeBC_with_SAB_RPL.csv")
pvals.to_csv(OUT_DIR / "correlation_p_beforeBC_with_SAB_RPL.csv")
df_use.to_csv(OUT_DIR / "input_used_beforeBC_with_SAB_RPL.csv", index=False)

print("\n✅ Saved BEFORE-BC correlation tables to:")
print(OUT_DIR)


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

# import pandas as pd
# import numpy as np
# from scipy.stats import pearsonr
# from pathlib import Path

# # -----------------------------
# # Paths
# # -----------------------------
# PCA_FILE   = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/PCA-after BC outputs/PCA_components_afterBC.csv")
# PHENO_FILE = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv")
# NAME_FILE  = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/4. Cluster/1. All-sig/Name Match.xlsx")
# DEMO_FILE  = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Demographics_EPL_final.xlsx")

# OUT_DIR = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/correlation output")
# OUT_DIR.mkdir(parents=True, exist_ok=True)

# # -----------------------------
# # Load files
# # -----------------------------
# pca = pd.read_csv(PCA_FILE)
# pheno = pd.read_csv(PHENO_FILE)
# name_map = pd.read_excel(NAME_FILE)
# demo = pd.read_excel(DEMO_FILE)

# # -----------------------------
# # Clean column names
# # -----------------------------
# for df in [pca, pheno, name_map, demo]:
#     df.columns = df.columns.str.replace("\u00A0", " ").str.strip()

# # Fix demographics typo
# if "RPL " in demo.columns and "RPL" not in demo.columns:
#     demo = demo.rename(columns={"RPL ": "RPL"})

# # Drop junk index column
# pca = pca.drop(columns=["Unnamed: 0"], errors="ignore")

# # -----------------------------
# # Attach Sample Name by row order
# # -----------------------------
# n = min(len(pca), len(pheno))
# pca = pca.reset_index(drop=True)
# pheno = pheno.reset_index(drop=True)

# pca.loc[:n-1, "Sample Name"] = pheno.loc[:n-1, "Sample Name"].astype(str)

# # -----------------------------
# # Map Sample Name -> StudyID
# # -----------------------------
# pca["SampleName_d"] = pca["Sample Name"].astype(str) + ".d"

# name_dict = dict(
#     zip(
#         name_map["Sample Name"].astype(str),
#         name_map["Sample Name_2"].astype(str)
#     )
# )

# pca["StudyID"] = pca["SampleName_d"].map(name_dict).fillna(pca["SampleName_d"])

# # -----------------------------
# # Merge demographics
# # -----------------------------
# pca = pca.merge(
#     demo[["StudyID", "SAB", "RPL"]],
#     on="StudyID",
#     how="left"
# )

# # -----------------------------
# # Define Sample Type (LOSS vs CONTROL)
# # -----------------------------
# # 1 = Loss (RPL), 0 = Control
# pca["Sample Type"] = pca["RPL"]

# print("\nSample Type counts BEFORE filtering:")
# print(pca["Sample Type"].value_counts(dropna=False))

# # -----------------------------
# # Encode Batch (optional)
# # -----------------------------
# def factorize_series(s):
#     return pd.factorize(s.astype(str))[0] + 1

# if "color_batch" in pca.columns:
#     pca["Batch"] = factorize_series(pca["color_batch"])

# # Ensure numeric
# for col in ["Sample Type", "Batch", "RPL", "SAB"]:
#     if col in pca.columns:
#         pca[col] = pd.to_numeric(pca[col], errors="coerce")

# # -----------------------------
# # Select variables for correlation
# # -----------------------------
# pc_cols = [c for c in pca.columns if c.startswith("PC")]

# df_use = pca[
#     pc_cols + ["Sample Type", "Batch", "RPL", "SAB"]
# ].copy()

# # -----------------------------
# # Rename PCs to PC1, PC2, ...
# # -----------------------------
# df_use = df_use.rename(
#     columns={f"PC{i}": f"PC{i+1}" for i in range(200)}
# )

# # Recompute PC column list AFTER renaming
# pc_cols = [c for c in df_use.columns if c.startswith("PC")]

# # -----------------------------
# # Drop rows missing PCs ONLY
# # -----------------------------
# df_use = df_use.dropna(subset=pc_cols)

# # -----------------------------
# # Drop samples without Sample Type (recommended)
# # -----------------------------
# df_use = df_use.dropna(subset=["Sample Type"])

# print("\nSample Type counts USED in correlation:")
# print(df_use["Sample Type"].value_counts())

# assert df_use["Sample Type"].nunique() > 1, \
#     "ERROR: Sample Type has no variation after filtering!"

# print("AFTER BC rows used:", len(df_use))

# # -----------------------------
# # Compute correlations
# # -----------------------------
# cols = df_use.columns.tolist()

# corr = pd.DataFrame(index=cols, columns=cols, dtype=float)
# pvals = pd.DataFrame(index=cols, columns=cols, dtype=float)

# for c1 in cols:
#     for c2 in cols:
#         r, p = pearsonr(df_use[c1], df_use[c2])
#         corr.loc[c1, c2] = r
#         pvals.loc[c1, c2] = p

# # -----------------------------
# # Save outputs
# # -----------------------------
# corr.to_csv(OUT_DIR / "correlation_r_afterBC_with_SAB_RPL.csv")
# pvals.to_csv(OUT_DIR / "correlation_p_afterBC_with_SAB_RPL.csv")
# df_use.to_csv(OUT_DIR / "input_used_afterBC_with_SAB_RPL.csv", index=False)

# print("\n✅ Saved AFTER-BC correlation tables to:")
# print(OUT_DIR)


