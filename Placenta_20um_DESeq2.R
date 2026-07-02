# =============================================================================
# Placenta 20uM PFBS — DESeq2 Differential Expression Analysis
# Skill references: deseq2-basics, de-results, de-visualization
# Key fix vs. prior Rmd: paired design (~ sample_id + treatment) blocks
# between-donor variance, which is the dominant noise source in human tissue.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
library(DESeq2)
library(here)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(purrr)
library(apeglm)

# ---- 1. LOAD COUNT DATA -----------------------------------------------------

base_25 <- here("data", "placenta_2025")
base_26 <- here("data", "placenta_2026")

load_tab <- function(path, col_name) {
  df <- read.csv(path, header = TRUE, sep = "\t")
  colnames(df)[2] <- col_name
  df
}

# 2025 donors (study = first)
n1_20um_25  <- load_tab(file.path(base_25, "1_20um_counts.tabular"),    "n1_20um_25")
n1_ctrl_25  <- load_tab(file.path(base_25, "1_Control_counts.tabular"), "n1_ctrl_25")
n2_20um_25  <- load_tab(file.path(base_25, "2_20um_counts.tabular"),    "n2_20um_25")
n2_ctrl_25  <- load_tab(file.path(base_25, "2_Control_counts.tabular"), "n2_ctrl_25")
n3_20um_25  <- load_tab(file.path(base_25, "3_20um_counts.tabular"),    "n3_20um_25")
n3_ctrl_25  <- load_tab(file.path(base_25, "3_Control_counts.tabular"), "n3_ctrl_25")

# 2026 donors (study = second)
n1_20um_26  <- load_tab(file.path(base_26, "n1_20um_26 featureCounts.tabular"),    "n1_20um_26")
n1_ctrl_26  <- load_tab(file.path(base_26, "n1_control_26 featureCounts.tabular"), "n1_ctrl_26")
n2_20um_26  <- load_tab(file.path(base_26, "n2_20um_26 featureCounts.tabular"),    "n2_20um_26")
n2_ctrl_26  <- load_tab(file.path(base_26, "n2_control_26 featureCounts.tabular"), "n2_ctrl_26")
n3_20um_26  <- load_tab(file.path(base_26, "n3_20um_26 featureCounts.tabular"),    "n3_20um_26")
n3_ctrl_26  <- load_tab(file.path(base_26, "n3_control_26 featureCounts.tabular"), "n3_ctrl_26")
n5_20um_26  <- load_tab(file.path(base_26, "n5_20um_26 featureCounts.tabular"),    "n5_20um_26")
n5_ctrl_26  <- load_tab(file.path(base_26, "n5_control_26 featureCounts.tabular"), "n5_ctrl_26")

# Merge all 14 samples
df_list <- list(
  n1_20um_25, n1_ctrl_25,
  n2_20um_25, n2_ctrl_25,
  n3_20um_25, n3_ctrl_25,
  n1_20um_26, n1_ctrl_26,
  n2_20um_26, n2_ctrl_26,
  n3_20um_26, n3_ctrl_26,
  n5_20um_26, n5_ctrl_26
)

count_data <- purrr::reduce(df_list, dplyr::full_join, by = "Geneid") %>%
  column_to_rownames(var = "Geneid")
count_data[is.na(count_data)] <- 0

cat("Count matrix dimensions:", nrow(count_data), "genes x", ncol(count_data), "samples\n")
stopifnot(ncol(count_data) == 14)

# ---- 2. METADATA ------------------------------------------------------------

meta_data <- data.frame(
  treatment = factor(c("treated","control","treated","control","treated","control",
                       "treated","control","treated","control","treated","control",
                       "treated","control")),
  sample_id = factor(c("one","one","two","two","three","three",
                       "four","four","five","five","six","six",
                       "eight","eight")),
  sex       = factor(c("M","M","M","M","M","M","F","F","F","F","M","M","F","F")),
  study     = factor(c("first","first","first","first","first","first",
                       "second","second","second","second","second","second",
                       "second","second")),
  row.names = colnames(count_data)
)

# Verify alignment
stopifnot(all(colnames(count_data) == rownames(meta_data)))
cat("Metadata/count matrix aligned: OK\n")
print(meta_data)

# =============================================================================
# SECTION 1: DATA DIAGNOSTICS (before any modeling)
# =============================================================================

out_fig <- here("output", "figures")
out_pl  <- here("output", "placenta")

# -- 1a. Library sizes --------------------------------------------------------
lib_sizes <- colSums(count_data)
lib_mean  <- mean(lib_sizes)
lib_sd    <- sd(lib_sizes)
outlier   <- lib_sizes > lib_mean + 2 * lib_sd | lib_sizes < lib_mean - 2 * lib_sd

pdf(file.path(out_fig, "diag_library_sizes.pdf"), width = 10, height = 5)
barplot(lib_sizes,
        las = 2, cex.names = 0.7,
        col = ifelse(outlier, "red", "steelblue"),
        main = "Library Sizes per Sample",
        ylab = "Total Counts")
abline(h = lib_mean, lty = 2, col = "red")
abline(h = c(lib_mean - 2*lib_sd, lib_mean + 2*lib_sd), lty = 3, col = "orange")
legend("topright", legend = c("Mean", "Mean ± 2SD", "Outlier"),
       lty = c(2,3,NA), col = c("red","orange","red"),
       pch = c(NA,NA,15), pt.cex = 1.5, bty = "n")
dev.off()

cat("\n--- Library Sizes ---\n")
print(sort(lib_sizes))
if (any(outlier)) {
  cat("WARNING: Potential outlier samples (>mean ± 2SD):", names(lib_sizes)[outlier], "\n")
} else {
  cat("Library sizes: no extreme outliers detected\n")
}

# -- 1b. Raw count distribution (log scale boxplot) ---------------------------
pdf(file.path(out_fig, "diag_count_distribution.pdf"), width = 12, height = 5)
log_counts <- log2(count_data + 1)
boxplot(log_counts, las = 2, cex.axis = 0.7,
        col = ifelse(meta_data$treatment == "treated", "#d95f02", "#1b9e77"),
        main = "log2(count+1) Distribution per Sample",
        ylab = "log2(count + 1)")
legend("topright", legend = c("Treated", "Control"),
       fill = c("#d95f02", "#1b9e77"), bty = "n")
dev.off()

# -- 1c. Zero proportion per gene ---------------------------------------------
zero_prop <- rowMeans(count_data == 0)
pct_half_zero <- round(mean(zero_prop > 0.5) * 100, 1)

pdf(file.path(out_fig, "diag_zero_proportions.pdf"), width = 7, height = 5)
hist(zero_prop, breaks = 50, col = "steelblue", border = "white",
     main = "Proportion of Zero Counts per Gene",
     xlab = "Proportion of samples with zero count")
abline(v = 0.5, col = "red", lty = 2)
text(0.52, par("usr")[4] * 0.9,
     paste0(pct_half_zero, "% genes >50% zeros"),
     col = "red", adj = 0, cex = 0.85)
dev.off()

cat("\nGenes with >50% zero counts:", sum(zero_prop > 0.5),
    paste0("(", pct_half_zero, "% of total)\n"))

# -- 1d. Mitochondrial gene check ---------------------------------------------
mt_genes <- grep("^MT-", rownames(count_data), value = TRUE)
cat("\nMitochondrial genes detected:", length(mt_genes), "\n")
if (length(mt_genes) > 0) {
  mt_counts <- colSums(count_data[mt_genes, , drop = FALSE])
  mt_pct    <- round(mt_counts / lib_sizes * 100, 1)
  cat("MT gene % of library per sample:\n")
  print(mt_pct)
  if (any(mt_pct > 20)) {
    cat("WARNING: MT% > 20% in some samples — may indicate low RNA quality\n")
  }
}

# Top 20 most expressed genes across all samples
top20 <- sort(rowMeans(count_data), decreasing = TRUE)[1:20]
cat("\nTop 20 most expressed genes (mean across samples):\n")
print(round(top20))

# =============================================================================
# SECTION 2: CREATE DESeqDataSet + PRE-FILTER
# Skill: deseq2-basics — paired design and group-aware filtering
# =============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = meta_data,
  design    = ~ sample_id + treatment   # paired design: key fix
)

# Group-aware pre-filter: require >= 10 counts in at least the smallest group
# (7 treated vs 7 control, so min group = 7)
min_group <- min(table(dds$treatment))
keep <- rowSums(counts(dds) >= 10) >= min_group
dds  <- dds[keep, ]
cat("\nGenes after group-aware filtering:", nrow(dds),
    "(removed", sum(!keep), "low-count genes)\n")

# Set reference levels
dds$treatment <- relevel(dds$treatment, ref = "control")
dds$sample_id <- relevel(dds$sample_id, ref = "one")

# =============================================================================
# SECTION 3: PRE-MODEL DIAGNOSTICS (VST, blind=TRUE)
# Skill: de-visualization — PCA and sample distance heatmap
# =============================================================================

vsd <- vst(dds, blind = TRUE)

# -- 3a. PCA plots ------------------------------------------------------------
pca_data   <- plotPCA(vsd, intgroup = c("treatment", "study"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pdf(file.path(out_fig, "diag_pca_premodel.pdf"), width = 10, height = 8)

# By treatment
p_treat <- ggplot(pca_data, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("control" = "#1b9e77", "treated" = "#d95f02")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Treatment") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# By study (batch)
p_study <- ggplot(pca_data, aes(PC1, PC2, color = study)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("first" = "#7570b3", "second" = "#e7298a")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Study (Batch)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# By sex
pca_sex <- plotPCA(vsd, intgroup = "sex", returnData = TRUE)
p_sex <- ggplot(pca_sex, aes(PC1, PC2, color = sex)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("M" = "#2166ac", "F" = "#d6604d")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Sex") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# By sample_id (donor)
pca_donor <- plotPCA(vsd, intgroup = "sample_id", returnData = TRUE)
p_donor <- ggplot(pca_donor, aes(PC1, PC2, color = sample_id)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Donor") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

print(p_treat)
print(p_study)
print(p_sex)
print(p_donor)
dev.off()

# -- 3b. Sample distance heatmap ----------------------------------------------
sampleDists <- dist(t(assay(vsd)))
distMatrix  <- as.matrix(sampleDists)

ann_colors <- list(
  treatment = c(control = "#1b9e77", treated = "#d95f02"),
  study     = c(first = "#7570b3", second = "#e7298a"),
  sex       = c(M = "#2166ac", F = "#d6604d")
)

pdf(file.path(out_fig, "diag_sample_distances.pdf"), width = 9, height = 8)
pheatmap(distMatrix,
         annotation_col = meta_data[, c("treatment", "study", "sex")],
         annotation_colors = ann_colors,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         main = "Sample Distance Matrix")
dev.off()

# -- 3c. Size factors ---------------------------------------------------------
dds <- estimateSizeFactors(dds)
cat("\nSize factors:\n")
print(round(sizeFactors(dds), 3))
sf_range <- max(sizeFactors(dds)) / min(sizeFactors(dds))
if (sf_range > 3) {
  cat("WARNING: Size factor range is", round(sf_range, 1),
      "— large variation may indicate library quality issues\n")
}

# =============================================================================
# SECTION 4: RUN DESeq2
# Skill: deseq2-basics — standard workflow, independent filtering ENABLED
# =============================================================================

cat("\nRunning DESeq2...\n")
dds <- DESeq(dds)

cat("\nResult names available:\n")
print(resultsNames(dds))

# Extract results — independent filtering ENABLED (do NOT pass independentFiltering=FALSE)
res <- results(dds, alpha = 0.05)

# LFC shrinkage for visualization and ranking (apeglm recommended by deseq2-basics skill)
resLFC <- lfcShrink(dds,
                    coef = "treatment_treated_vs_control",
                    type = "apeglm")

# =============================================================================
# SECTION 5: POST-MODEL DIAGNOSTICS
# Skill: de-visualization — dispersion, p-value histogram, MA plot
# =============================================================================

pdf(file.path(out_fig, "diag_postmodel.pdf"), width = 10, height = 12)
par(mfrow = c(2, 2))

# Dispersion plot
plotDispEsts(dds, main = "Dispersion Estimates")

# P-value histogram (un-shrunken — use res not resLFC)
res_df <- as.data.frame(res)
hist(res_df$pvalue, breaks = 50, col = "steelblue", border = "white",
     main = "P-value Distribution\n(uniform + spike at 0 = good)",
     xlab = "p-value")

# MA plot with shrunken LFCs
plotMA(resLFC, ylim = c(-4, 4),
       main = "MA Plot (shrunken LFC)\nblue = padj < 0.05")

# Cook's distances (outlier check)
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2, cex.axis = 0.6,
        main = "Cook's Distances per Sample",
        ylab = "log10(Cook's distance)")

dev.off()

# =============================================================================
# SECTION 6: RESULTS
# Skill: de-results — summary, filtering, export
# =============================================================================

cat("\n========== DESeq2 Results Summary ==========\n")
summary(res, alpha = 0.05)

# Add gene column and reorder
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene","log2FoldChange","baseMean","lfcSE","stat","pvalue","padj")]

# Significant genes: padj < 0.05
sig_padj <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

# Stringent: padj < 0.05 AND |LFC| > 1
sig_strict <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

cat("\nSignificant genes (padj < 0.05):               ", nrow(sig_padj), "\n")
cat("Significant genes (padj < 0.05 & |LFC| > 1):  ", nrow(sig_strict), "\n")

if (nrow(sig_padj) > 0) {
  cat("\nTop 10 significant genes:\n")
  print(head(sig_padj, 10))
}

# -- Export -------------------------------------------------------------------
write.csv(res_df,
          file = file.path(out_pl, "differential_expression_20um_paired.csv"),
          row.names = FALSE)

if (nrow(sig_padj) > 0) {
  write.csv(sig_padj,
            file = file.path(out_pl, "significant_genes_20um_paired.csv"),
            row.names = FALSE)
}

write.csv(as.data.frame(count_data),
          file = file.path(out_pl, "combined_count_20um.csv"),
          row.names = TRUE)

cat("\nOutputs saved to:", out_pl, "\n")
cat("Diagnostic plots saved to:", out_fig, "\n")
cat("\nDone.\n")
