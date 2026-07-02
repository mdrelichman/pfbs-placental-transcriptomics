# =============================================================================
# Over-Representation Analysis: HIF1A Targets in Female Placenta 5 uM PFBS DGE
# Method: Hypergeometric test via phyper()
# Gene set: Harmonizome ChEA HIF1A Transcription Factor Targets (314 genes)
#           fetched live from the Harmonizome REST API
# =============================================================================

library(tidyverse)
library(here)
library(jsonlite)

# ---- 1. LOAD DGE RESULTS ----------------------------------------------------

dge <- read.csv(
  here("output", "placenta", "female_differential_expression_5um_26.csv"),
  row.names = 1
)

# Universe = all genes tested in this DGE
universe <- dge$gene

# Biologically significant genes: padj < 0.05 AND |log2FC| >= 1
sig_genes <- dge %>%
  filter(pvalue < 0.05) %>%
  pull(gene)

cat("Universe size:", length(universe), "\n")
cat("Significant genes (pvalue < 0.05):", length(sig_genes), "\n")

# ---- 2. FETCH HARMONIZOME HIF1A GENE SET VIA API ----------------------------
# Pulls all 314 target genes directly — no manual download needed.
# API docs: https://maayanlab.cloud/Harmonizome/api/1.0

api_url <- paste0(
  "https://maayanlab.cloud/Harmonizome/api/1.0/",
  "gene_set/HIF1A/ChEA+Transcription+Factor+Targets"
)

cat("Fetching HIF1A gene set from Harmonizome API...\n")
response       <- fromJSON(api_url, flatten = TRUE)
hif1a_genes_raw <- response$associations$gene.symbol

cat("Genes retrieved from API:", length(hif1a_genes_raw), "\n")

# Restrict gene set to genes present in the tested universe
hif1a_in_universe <- intersect(hif1a_genes_raw, universe)

cat("HIF1A Harmonizome targets (raw):", length(hif1a_genes_raw), "\n")
cat("HIF1A targets present in universe:", length(hif1a_in_universe), "\n")

# ---- 3. COMPUTE OVERLAP -----------------------------------------------------

overlap <- intersect(sig_genes, hif1a_in_universe)
cat("Overlap (sig genes ∩ HIF1A targets):", length(overlap), "\n")
cat("Overlapping genes:", paste(overlap, collapse = ", "), "\n")

# ---- 3b. DIRECTION OF OVERLAP GENES -----------------------------------------

overlap_direction <- dge %>%
  filter(gene %in% overlap) %>%
  select(gene, log2FoldChange, pvalue, padj) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "upregulated", "downregulated")) %>%
  arrange(log2FoldChange)

cat("\nDirection of overlapping HIF1A target genes:\n")
print(overlap_direction)

cat("\nSummary:\n")
print(table(overlap_direction$direction))

# ---- 4. HYPERGEOMETRIC TEST (phyper) ----------------------------------------
# phyper(q, m, n, k, lower.tail = FALSE)
#   q = overlap size - 1      (observed successes, minus 1 for upper tail)
#   m = HIF1A targets in universe  (white balls in urn)
#   n = universe - m               (black balls in urn)
#   k = significant genes          (balls drawn)

q <- length(overlap) - 1
m <- length(hif1a_in_universe)
n <- length(universe) - m
k <- length(sig_genes)

p_value <- phyper(q, m, n, k, lower.tail = FALSE)

# Fold enrichment: observed overlap / expected overlap
expected_overlap <- k * (m / length(universe))
fold_enrichment  <- length(overlap) / expected_overlap

cat("\n--- ORA Results ---\n")
cat("q (overlap - 1):", q, "\n")
cat("m (HIF1A targets in universe):", m, "\n")
cat("n (non-HIF1A genes in universe):", n, "\n")
cat("k (significant genes drawn):", k, "\n")
cat("Observed overlap:", length(overlap), "\n")
cat("Expected overlap:", round(expected_overlap, 3), "\n")
cat("Fold enrichment:", round(fold_enrichment, 3), "\n")
cat("Hypergeometric p-value:", p_value, "\n")

# ---- 5. SUMMARY TABLE -------------------------------------------------------

overlap_gene_names <- paste(overlap, collapse = "; ")

results_table <- tibble(
  gene_set          = "Harmonizome HIF1A TF Targets",
  universe_size     = length(universe),
  gene_set_in_univ  = m,
  sig_genes         = k,
  overlap           = length(overlap),
  expected_overlap  = round(expected_overlap, 3),
  fold_enrichment   = round(fold_enrichment, 3),
  p_value           = p_value,
  overlap_genes     = overlap_gene_names
)

print(results_table)

write.csv(
  results_table,
  here("output", "placenta", "ORA_HIF1A_female_5um.csv"),
  row.names = FALSE
)
