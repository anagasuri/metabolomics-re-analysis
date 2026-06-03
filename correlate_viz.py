import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------------
# Output directory (BEFORE BC)
# -----------------------------
OUT_DIR = Path(
    "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/Before BC/correlation output"
)

# -----------------------------
# Load correlation matrices
# -----------------------------
r_mat = pd.read_csv(
    OUT_DIR / "correlation_r_beforeBC_with_SAB_RPL.csv",
    index_col=0
)

p_mat = pd.read_csv(
    OUT_DIR / "correlation_p_beforeBC_with_SAB_RPL.csv",
    index_col=0
)

# -----------------------------
# Subset variables
# -----------------------------
keep = ["PC1", "PC2", "PC3", "Sample Type", "Batch", "RPL", "SAB"]

r = r_mat.loc[keep, keep]
p = p_mat.loc[keep, keep]

# -----------------------------
# Mask upper triangle
# -----------------------------
mask = np.triu(np.ones_like(r, dtype=bool))
sns.set_style("white")

# -----------------------------
# Significance annotations (P ONLY)
# -----------------------------
sig_annot = pd.DataFrame("", index=r.index, columns=r.columns)

for i in r.index:
    for j in r.columns:
        if p.loc[i, j] < 0.01:
            sig_annot.loc[i, j] = "*"

# =============================
# R-VALUE HEATMAP
# =============================
plt.figure(figsize=(7, 6))
sns.heatmap(
    r,
    mask=mask,
    vmin=-1, vmax=1,
    cmap="coolwarm",
    cbar_kws={"label": "Correlation Coefficient"},
    square=True
)

plt.title(
    "Correlation Matrix (R-values)\n |R| > 0.5",
    fontsize=14,
    pad=14
)
plt.xticks(rotation=0)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(
    OUT_DIR / "PC1_PC3_SampleType_Batch_RPL_SAB_R_sigP_beforeBC.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# =============================
# P-VALUE HEATMAP
# =============================
plt.figure(figsize=(7, 6))
sns.heatmap(
    p,
    mask=mask,
    annot=sig_annot,
    vmin=0, vmax=0.05,
    cmap="Blues_r",
    fmt="",
    annot_kws={"size": 18, "weight": "bold"},
    cbar_kws={"label": "P-value"},
    square=True
)
plt.title(
    "Correlation Matrix (P-values)\n* = p < 0.01",
    fontsize=14,
    pad=14
)
plt.xticks(rotation=0)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(
    OUT_DIR / "PC1_PC3_SampleType_Batch_RPL_SAB_P_sigP_beforeBC.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Saved BEFORE-BC correlation plots (stars = p < 0.01) to:", OUT_DIR)


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

# import pandas as pd
# import seaborn as sns
# import numpy as np
# import matplotlib.pyplot as plt
# from pathlib import Path

# # -----------------------------
# # Output directory (AFTER BC)
# # -----------------------------
# OUT_DIR = Path(
#     "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/PCA Batch Removal - REDO/correlation output"
# )

# # -----------------------------
# # Load correlation matrices
# # -----------------------------
# r_mat = pd.read_csv(
#     OUT_DIR / "correlation_r_afterBC_with_SAB_RPL.csv",
#     index_col=0
# )

# p_mat = pd.read_csv(
#     OUT_DIR / "correlation_p_afterBC_with_SAB_RPL.csv",
#     index_col=0
# )

# # -----------------------------
# # Variables to include
# # -----------------------------
# keep = ["PC1", "PC2", "PC3", "Sample Type", "Batch", "RPL", "SAB"]

# r = r_mat.loc[keep, keep]
# p = p_mat.loc[keep, keep]

# # -----------------------------
# # Mask upper triangle
# # -----------------------------
# mask = np.triu(np.ones_like(r, dtype=bool))
# sns.set_style("white")

# # -----------------------------
# # Significance annotations (P ONLY)
# # -----------------------------
# sig_annot = pd.DataFrame("", index=r.index, columns=r.columns)

# for i in r.index:
#     for j in r.columns:
#         if p.loc[i, j] < 0.01:
#             sig_annot.loc[i, j] = "*"

# # =============================
# # R-VALUE HEATMAP
# # =============================
# plt.figure(figsize=(7, 6))
# sns.heatmap(
#     r,
#     mask=mask,
#     vmin=-1, vmax=1,
#     cmap="coolwarm",
#     cbar_kws={"label": "Correlation Coefficient"},
#     square=True
# )

# plt.title(
#     "Correlation Matrix (R-values)\n |R| > 0.5",
#     fontsize=14,
#     pad=14
# )
# plt.xticks(rotation=0)
# plt.yticks(rotation=0)
# plt.tight_layout()
# plt.savefig(
#     OUT_DIR / "PC1_PC3_SampleType_Batch_RPL_SAB_R_sigP_afterBC.png",
#     dpi=300,
#     bbox_inches="tight"
# )
# plt.close()

# # =============================
# # P-VALUE HEATMAP
# # =============================
# plt.figure(figsize=(7, 6))
# sns.heatmap(
#     p,
#     mask=mask,
#     annot=sig_annot,
#     vmin=0, vmax=0.05,
#     cmap="Blues_r",
#     fmt="",
#     annot_kws={"size": 18, "weight": "bold"},
#     cbar_kws={"label": "P-value"},
#     square=True
# )
# plt.title(
#     "Correlation Matrix (P-values)\n* = p < 0.01",
#     fontsize=14,
#     pad=14
# )
# plt.xticks(rotation=0)
# plt.yticks(rotation=0)
# plt.tight_layout()
# plt.savefig(
#     OUT_DIR / "PC1_PC3_SampleType_Batch_RPL_SAB_P_sigP_afterBC.png",
#     dpi=300,
#     bbox_inches="tight"
# )
# plt.close()

# print("Saved AFTER-BC correlation plots (stars = p < 0.01) to:", OUT_DIR)

