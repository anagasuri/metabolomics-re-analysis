# -*- coding: utf-8 -*-
"""
Before BC PCA — replicate averaging + boxplots
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import seaborn as sns
from statannotations.Annotator import Annotator
from matplotlib.lines import Line2D

# ----------------------------------------------------------
# OUTPUT DIRECTORY
# ----------------------------------------------------------
OUT = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/PCA-BC-AVERAGED outputs")
OUT.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------
# LOAD FILES
# ----------------------------------------------------------
df1 = pd.read_csv(
    "/Users/amritanagasuri/Downloads/Q-TOF Data Analysis/2. imput-PCA/Before BC/neg,pos,pheno_inputs/Neg_Clean_NEWWITHREPS.csv"
)
df2 = pd.read_csv(
    "/Users/amritanagasuri/Downloads/Q-TOF Data Analysis/2. imput-PCA/Before BC/neg,pos,pheno_inputs/Pos_Clean_NEWWITHREPS.csv"
)
df3 = pd.read_csv(
    "/Users/amritanagasuri/Downloads/Q-TOF Data Analysis/2. imput-PCA/Before BC/neg,pos,pheno_inputs/pheno-clean-NEWWITHREPS.csv"
)

print("Loaded Neg:", df1.shape, "Pos:", df2.shape, "Pheno:", df3.shape)

# ----------------------------------------------------------
# REPLICATE AVERAGING
# ----------------------------------------------------------
def average_replicates(df):
    sample_cols = [c for c in df.columns if c.startswith("BI")]
    groups = {}

    for col in sample_cols:
        base = col.replace("_r1", "").replace("_r2", "")
        groups.setdefault(base, []).append(col)

    averaged = {base: df[cols].mean(axis=1) for base, cols in groups.items()}
    meta_cols = [c for c in df.columns if c not in sample_cols]

    return pd.concat([df[meta_cols], pd.DataFrame(averaged)], axis=1)

df1_avg = average_replicates(df1)
df2_avg = average_replicates(df2)
df = pd.concat([df1_avg, df2_avg], axis=0)

print("After replicate averaging + concat:", df.shape)

# ----------------------------------------------------------
# BUILD BATCH + PHENO MAPPINGS
# ----------------------------------------------------------
df3_clean = df3.copy()
df3_clean["base"] = df3_clean["Sample Name"].str.replace(r"_r[12]$", "", regex=True)

batch_map = df3_clean.set_index("base")["batch"].to_dict()
age_map   = df3_clean.set_index("base")["age"].to_dict()

column_mapping = {
    col: f"{batch_map[col]}_{age_map[col]}_{col}"
    for col in df.columns
    if col in batch_map
}

df.rename(columns=column_mapping, inplace=True)

# ----------------------------------------------------------
# PCA INPUT MATRIX
# ----------------------------------------------------------
df = df.set_index("chem_id")

loc1 = df.loc[:, df.columns.str.match(r"^1_")]
loc2 = df.loc[:, df.columns.str.match(r"^2_")]
loc3 = df.loc[:, df.columns.str.match(r"^3_")]

dr = pd.concat([loc1, loc2, loc3], axis=1).T

# Add color labels
dr["color_MC"] = np.where(dr.index.str.contains("_P"), "red", "black")
dr["color_batch"] = np.where(
    dr.index.str.contains("1_"),
    "#FFD700",
    np.where(dr.index.str.contains("2_"), "#2E8B57", "#006400"),
)

dr12 = dr[["color_MC", "color_batch"]]
dr11 = dr.drop(["color_MC", "color_batch"], axis=1)

# ----------------------------------------------------------
# CLEAN NANS + PCA
# ----------------------------------------------------------
dr11_cleaned = dr11.dropna(axis=1)
X_std = StandardScaler().fit_transform(dr11_cleaned)

pca = PCA(n_components=20)
pc = pca.fit_transform(X_std)

PCA_components = pd.DataFrame(pc, columns=[f"PC{i+1}" for i in range(pc.shape[1])])
PCA_components = pd.concat([PCA_components, dr12.reset_index(drop=True)], axis=1)

# ----------------------------------------------------------
# SHARED LEGEND HANDLES
# ----------------------------------------------------------
legend_MC = [
    Line2D([0], [0], marker='o', color='w', label='Loss',
           markerfacecolor='red', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Control',
           markerfacecolor='black', markersize=8)
]

legend_batch = [
    Line2D([0], [0], marker='o', color='w', label='1',
           markerfacecolor='#FFD700', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='2',
           markerfacecolor='#2E8B57', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='3',
           markerfacecolor='#006400', markersize=8)
]

# ----------------------------------------------------------
# PLOTS
# ----------------------------------------------------------

# === PCA variance explained (NO LEGEND) ===
plt.figure(figsize=(8, 5))
pcs = np.arange(1, pca.n_components_ + 1)
plt.bar(pcs, pca.explained_variance_ratio_)
plt.xticks(pcs, [f"PC{i}" for i in pcs], rotation=45)
plt.xlabel("Principal Components")
plt.ylabel("Variance Explained")
plt.savefig(OUT / "PCA_variance_bar.png", dpi=400, bbox_inches="tight")
plt.close()

# === PC1 vs PC2 (sample type) ===
plt.figure()
plt.scatter(PCA_components["PC1"], PCA_components["PC2"],
            c=PCA_components["color_MC"])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(handles=legend_MC, title="Groups",
           loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig(OUT / "PC1_vs_PC2_sample_type.png", dpi=400, bbox_inches="tight")
plt.close()

# === PC1 vs PC2 (batch) ===
plt.figure()
plt.scatter(PCA_components["PC1"], PCA_components["PC2"],
            c=PCA_components["color_batch"])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(handles=legend_batch, title="Batch",
           loc="center left", bbox_to_anchor=(1, 0.5))
plt.savefig(OUT / "PC1_vs_PC2_batch.png", dpi=400, bbox_inches="tight")
plt.close()

# ----------------------------------------------------------
# BOXPLOTS
# ----------------------------------------------------------

# === Boxplot PC1 by Batch ===
df_box = PCA_components.copy()
df_box["batch"] = df_box["color_batch"].replace(
    {"#FFD700": "1", "#2E8B57": "2", "#006400": "3"}
)

plt.figure(figsize=(8, 6))
ax = sns.boxplot(
    data=df_box,
    x="batch",
    y="PC1",
    order=["1", "2", "3"],
    palette={"1": "#FFD700", "2": "#2E8B57", "3": "#006400"},
)

pairs = [("1", "2"), ("2", "3"), ("1", "3")]
annotator = Annotator(ax, pairs, data=df_box, x="batch", y="PC1",
                      order=["1", "2", "3"])
annotator.configure(test="Mann-Whitney", text_format="star", loc="outside")
annotator.apply_and_annotate()

ax.legend(handles=legend_batch, title="Batch",
          loc="center left", bbox_to_anchor=(1, 0.5))

plt.xlabel("Batch")
plt.ylabel("PC1")
plt.savefig(OUT / "boxplot_PC1_by_batch.png", dpi=400, bbox_inches="tight")
plt.close()

# === Boxplot PC1 by Sample Type ===
df_box2 = PCA_components.copy()
df_box2["MC"] = df_box2["color_MC"].replace({"red": "Loss", "black": "Control"})

plt.figure(figsize=(8, 6))
ax = sns.boxplot(
    data=df_box2,
    x="MC",
    y="PC1",
    order=["Loss", "Control"],
    palette={"Loss": "red", "Control": "black"},
    medianprops={"color": "white", "linewidth": 1},
)

pairs2 = [("Loss", "Control")]
annotator = Annotator(ax, pairs2, data=df_box2,
                      x="MC", y="PC1", order=["Loss", "Control"])
annotator.configure(test="Mann-Whitney", text_format="star", loc="outside")
annotator.apply_and_annotate()

ax.legend(handles=legend_MC, title="Groups",
          loc="center left", bbox_to_anchor=(1, 0.5))

plt.xlabel("Sample Type")
plt.ylabel("PC1")
plt.savefig(OUT / "boxplot_PC1_by_sample_type.png", dpi=400, bbox_inches="tight")
plt.close()

print("🎉 Done — replicate-averaged PCA + boxplots saved to:", OUT)
