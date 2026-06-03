library(limma)
library(readr)
library(ggplot2)
library(ggpubr)

# -------------------------------------------------------------------------
# Define output directory (no working directory change)
# -------------------------------------------------------------------------
out_dir <- "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/2. imput-PCA/removeBatchEffect_limma"

# -------------------------------------------------------------------------
# 1. NEGATIVE MODE
# -------------------------------------------------------------------------
df1 <- read_csv("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/Neg_Clean_NEWWITHREPS.csv")
df2 <- read_csv("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv")

head(df1)
head(df2)

# Rename first column to "name"
colnames(df1)[1] <- "name"
row_names <- df1$name
df1 <- df1[, -1, drop = FALSE]
rownames(df1) <- row_names

# Align samples
common_samples <- intersect(colnames(df1), df2$`Sample Name`)
if (length(common_samples) == 0) stop("No common samples found between feature table and metadata")
df1 <- df1[, common_samples, drop = FALSE]
df2 <- df2[match(common_samples, df2$`Sample Name`), ]
stopifnot(all(colnames(df1) == df2$`Sample Name`))

# Transform + batch correction
X_log <- log1p(as.matrix(df1))
batch <- factor(df2$batch)
X_corrected <- removeBatchEffect(X_log, batch = batch)
X_final <- expm1(X_corrected)
rownames(X_final) <- row_names

# Save corrected output
corrected_df <- as.data.frame(X_final)
corrected_df <- cbind(name = rownames(corrected_df), corrected_df)
cat("\nFinal NEG data structure:\n")
str(corrected_df[1:5, 1:5])
write_csv(corrected_df, file.path(out_dir, "Neg_removebatcheffect_corrected.csv"))

# -------------------------------------------------------------------------
# 2. POSITIVE MODE
# -------------------------------------------------------------------------
df3 <- read_csv("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/Pos_Clean_NEWWITHREPS.csv")
df4 <- read_csv("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/pheno-clean-NEWWITHREPS.csv")

colnames(df3)[1] <- "name"
row_names <- df3$name
df3 <- df3[, -1, drop = FALSE]
rownames(df3) <- row_names

common_samples <- intersect(colnames(df3), df4$`Sample Name`)
if (length(common_samples) == 0) stop("No common samples found between feature table and metadata")
df3 <- df3[, common_samples, drop = FALSE]
df4 <- df4[match(common_samples, df4$`Sample Name`), ]
stopifnot(all(colnames(df3) == df4$`Sample Name`))

X_log <- log1p(as.matrix(df3))
batch <- factor(df4$batch)
X_corrected <- removeBatchEffect(X_log, batch = batch)
X_final <- expm1(X_corrected)
rownames(X_final) <- row_names

corrected_df <- as.data.frame(X_final)
corrected_df <- cbind(name = rownames(corrected_df), corrected_df)
cat("\nFinal POS data structure:\n")
str(corrected_df[1:5, 1:5])
write_csv(corrected_df, file.path(out_dir, "Pos_removebatcheffect_corrected.csv"))

cat("\n✅ Saved corrected CSVs to:\n", out_dir, "\n")

# -------------------------------------------------------------------------
# 3. PCA + BOXPLOTS
# -------------------------------------------------------------------------
df3 <- read_csv("/Users/amritanagasuri/Desktop/new Q-TOF Analysis/final pos, neg, pheno/Pos_Clean_NEWWITHREPS.csv")
colnames(df3)[1] <- "name"
feature_matrix <- as.data.frame(df3[, -1])
rownames(feature_matrix) <- df3$name
X_original <- feature_matrix

df2 <- df4
i <- which(names(df2) == "Sample Name")
if (length(i) == 1) names(df2)[i] <- "sample_name"

# --- PCA plot function ---
plot_pca_only <- function(data_matrix, metadata, title) {
  common_samples <- intersect(colnames(data_matrix), metadata$sample_name)
  if (length(common_samples) == 0) stop("No common samples found between data matrix and metadata")
  data_matrix <- data_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$sample_name), ]
  stopifnot(all(colnames(data_matrix) == metadata$sample_name))
  pca_data <- prcomp(t(data_matrix), center = TRUE, scale. = TRUE)
  var_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2)
  pca_df <- as.data.frame(pca_data$x[, 1:2])
  pca_df$Batch <- factor(metadata$batch)
  pca_df$Sample <- rownames(pca_df)
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3) +
    labs(
      x = sprintf("PC1 (%.1f%%)", var_explained[1] * 100),
      y = sprintf("PC2 (%.1f%%)", var_explained[2] * 100),
      title = title
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.title.x = element_text(face = "plain"),
      axis.title.y = element_text(face = "plain"),
      axis.text = element_text(face = "plain"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.line = element_line(color = "black")
    )
}

# --- Boxplot function ---
plot_pc1_boxplot <- function(data_matrix, metadata) {
  common_samples <- intersect(colnames(data_matrix), metadata$sample_name)
  if (length(common_samples) == 0) stop("No common samples found between data matrix and metadata")
  data_matrix <- data_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, metadata$sample_name), ]
  pca_data <- prcomp(t(data_matrix), center = TRUE, scale. = TRUE)
  pc1_df <- data.frame(PC1 = pca_data$x[, 1], Batch = factor(metadata$batch))
  ggplot(pc1_df, aes(x = Batch, y = PC1, fill = Batch)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    labs(x = "Batch", y = "PC1 Score") +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = combn(levels(pc1_df$Batch), 2, simplify = FALSE),
      method.args = list(alternative = "two.sided", exact = FALSE),
      label = "p.signif"
    ) +
    scale_fill_brewer(palette = "Set1") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(face = "plain"),
      axis.title.y = element_text(face = "plain"),
      axis.text = element_text(face = "plain"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "plain"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.line = element_line(color = "black")
    )
}

# --- Generate plots ---
pca_original <- plot_pca_only(X_original, df2, "PCA Before Normalization")
pca_before   <- plot_pca_only(X_log, df2, "PCA After LOESS")
pca_after    <- plot_pca_only(X_corrected, df2, "PCA After Batch Correction")

box_original <- plot_pc1_boxplot(X_original, df2)
box_before   <- plot_pc1_boxplot(X_log, df2)
box_after    <- plot_pc1_boxplot(X_corrected, df2)

# --- Save plots to out_dir ---
ggsave(file.path(out_dir, "PCA_original.png"), pca_original, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "PCA_LOESS.png"), pca_before, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "PCA_after.png"), pca_after, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "PC1_box_original.png"), box_original, width = 6, height = 6, dpi = 300)
ggsave(file.path(out_dir, "PC1_box_LOESS.png"), box_before, width = 6, height = 6, dpi = 300)
ggsave(file.path(out_dir, "PC1_box_after.png"), box_after, width = 6, height = 6, dpi = 300)

cat("\n✅ All plots and CSVs saved to:\n", out_dir, "\n")

# --- Print summaries ---
print(pca_original)
print(box_original)
print(pca_before)
print(box_before)
print(pca_after)
print(box_after)
