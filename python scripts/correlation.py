# -*- coding: utf-8 -*-
"""
Before Batch Correction — Correlation Plots (FIXED)
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pathlib import Path 

OUT_DIR = Path("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/new plots correlation")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv('/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/PCA_components_beforeBC.csv')

df = df.drop(columns=['Unnamed: 0'], errors='ignore')

# convert PC columns to numeric (CRITICAL FIX)
pc_cols = [c for c in df.columns if c.startswith("PC")]
df[pc_cols] = df[pc_cols].apply(pd.to_numeric, errors="coerce")

# keep only PC columns (remove batch & sample type)
df = df[pc_cols]

# -----------------------------
# Correlation / p-values
# -----------------------------
correlation_matrix = df.corr(method='pearson')

p_values = pd.DataFrame(np.zeros((df.shape[1], df.shape[1])),
                        columns=df.columns, index=df.columns)

for col1 in df:
    for col2 in df:
        r, p = pearsonr(df[col1], df[col2])
        correlation_matrix.loc[col1, col2] = r
        p_values.loc[col1, col2] = p

mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))

sns.set_style("white")

# annotation
correlation_annotations = correlation_matrix.applymap(
    lambda x: f"{x:.2f}" if abs(x) > 0 else ""
)

p_value_annotations = p_values.applymap(
    lambda x: '*' if x < 0.05 and x >= 0.01 else '**' if x < 0.01 else ""
)

# -----------------------------
# Plot R-values
# -----------------------------
plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, mask=mask, annot=correlation_annotations,
            vmin=-1, vmax=1, cmap='coolwarm', fmt="",
            annot_kws={'size':12}, cbar_kws={'label':'Correlation Coefficient'})
plt.title('Correlation Matrix (R-values)', fontsize=18)
plt.savefig(OUT_DIR / "correlation_matrix.png", dpi=300, bbox_inches='tight')
plt.close()

# -----------------------------
# Plot P-values
# -----------------------------
plt.figure(figsize=(10, 8))
sns.heatmap(p_values, mask=mask, annot=p_value_annotations,
            vmin=0, vmax=0.05, cmap='Blues', fmt="",
            annot_kws={'size':12}, cbar_kws={'label':'P-value'})
plt.title('Correlation Matrix (P-values)', fontsize=18)
plt.savefig(OUT_DIR / "p_values_matrix.png", dpi=300, bbox_inches='tight')
plt.close()

# save CSV outputs
correlation_matrix.to_csv(OUT_DIR / "correlation_matrix.csv")
p_values.to_csv(OUT_DIR / "p_values_matrix.csv")
df.to_csv(OUT_DIR / "input_dataframe_used.csv")

print("All outputs saved to:", OUT_DIR)
