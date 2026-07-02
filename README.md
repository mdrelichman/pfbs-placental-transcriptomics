# PFBS Placental Transcriptomics — Thesis Analysis Scripts

Analysis code for a master's thesis investigating the transcriptomic effects of
**PFBS (perfluorobutane sulfonate)** exposure on the human placenta, using bulk
RNA-seq of both an *in vitro* trophoblast model (HTR-8/SVneo cells) and *ex vivo*
term placental tissue. The workflow covers transcript quantification, differential
expression (DESeq2), surrogate variable analysis, gene set enrichment (GSEA),
over-representation analysis (ORA), and cross-comparison visualization.

> **Note:** This repository contains **analysis scripts only.** Raw and processed
> sequencing data are not included here (see *Data availability* below). The scripts
> are published for transparency and reproducibility alongside the thesis/manuscript.

## Study design at a glance

- **Systems:** HTR-8/SVneo trophoblast cells (*in vitro*) and term placental tissue (*ex vivo*)
- **Exposures:** PFBS at 5 µM and 20 µM vs. vehicle control
- **Stratification:** placental tissue analyzed with attention to fetal sex
- **Question:** Which genes and biological pathways respond to PFBS exposure in placental cells and tissue?

## Repository contents

Scripts are grouped by role. Most are R Markdown (`.Rmd`); a few are plain R (`.R`).

### Differential expression (DESeq2)
| File | Purpose |
| --- | --- |
| `cells_PFBS.Rmd` | DE analysis of HTR-8/SVneo trophoblast cells, PFBS 5 µM vs. control. |
| `Placenta_20um_DESeq2.R` | DE analysis of term placental tissue, 20 µM (combined cohorts; paired design blocking between-donor variance, with QC diagnostics). |
| `Term_Placenta_Combined_5.Rmd` | Sex-stratified DE analysis of term placental tissue, 5 µM (combined cohorts, male/female modeled separately), **plus log2FC correlation analyses**: 5 µM vs 20 µM dose correlation within each sex, and HTR-8/SVneo cells vs male/female term tissue at 5 µM (`cor.test` + smoothScatter plots). |

### Gene set enrichment & over-representation
| File | Purpose |
| --- | --- |
| `Gene Set Enrichment Analysis (cells).Rmd` | GSEA of HTR-8/SVneo results (Hallmark, GO:BP, KEGG). |
| `Gene Set Enrichment Analysis (Placenta_5).Rmd` | GSEA of placental tissue, 5 µM. |
| `Gene Set Enrichment Analysis (Placenta_20).Rmd` | GSEA of placental tissue, 20 µM. |
| `GSEA_Sex_Comparison_Placenta.Rmd` | GSEA comparison stratified by fetal sex. |
| `ORA_HIF1A_female_5um.Rmd` | Over-representation analysis focused on HIF1A targets (female, 5 µM). |

### Comparison & visualization
| File | Purpose |
| --- | --- |
| `Comparisons and Correlations and Heatmap.Rmd` | Cross-condition correlations and expression heatmaps. |
| `GeneComparison.Rmd` | Overlap of differentially expressed genes across comparisons (UpSet / Venn). |

## Software requirements

- **R** (≥ 4.2 recommended) and **Bioconductor** (≥ 3.16).

**Bioconductor packages:** `DESeq2`, `apeglm`, `AnnotationDbi`, `org.Hs.eg.db`

**CRAN packages:** `here`, `tidyverse` (dplyr, ggplot2, stringr, tidyr, forcats, readr, purrr, tibble),
`readxl`, `openxlsx`, `pheatmap`, `RColorBrewer`, `ggrepel`, `UpSetR`, `ggVennDiagram`,
`DT`, `jsonlite`

Install example:

```r
install.packages(c("tidyverse","here","readxl","openxlsx","pheatmap","RColorBrewer",
                   "ggrepel","UpSetR","ggVennDiagram","DT","jsonlite"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","apeglm","AnnotationDbi","org.Hs.eg.db"))
```

## Running the scripts

1. Clone this repository and open **`pfbs-placental-transcriptomics.Rproj`** in RStudio.
   Opening the project sets the working directory to the repository root.
2. Restore the exact package versions used for the analysis (see *Reproducibility* below).
3. Recreate the expected data layout and place the input data:
   ```
   data/     <- input count matrices, gene sets, etc.
   output/   <- created by the scripts (results, figures); subfolders: output/cells, output/figures, output/placenta
   ```
4. Open a script and run it chunk-by-chunk, or knit the `.Rmd` to HTML.

### Reproducibility

- **Portable paths:** Scripts use the [`here`](https://here.r-lib.org/) package instead of
  absolute machine paths, e.g. `here("data", "cells", "N1_5um_counts.tabular")`. Because
  `here()` anchors to the project root, the same code runs on any computer that opens the
  `.Rproj`. *(Path conversion is being rolled out script-by-script; `cells_PFBS.Rmd` is the
  reference example.)*
- **Package versions:** This project uses [`renv`](https://rstudio.github.io/renv/). After
  opening the project, run `renv::restore()` to install the exact package versions recorded
  in `renv.lock`.
- **R version:** Developed with R ≥ 4.2. Run `sessionInfo()` for the full environment.

## Data availability

Raw and processed sequencing data are **not** stored in this repository. They are
available from the corresponding author on reasonable request / via the associated
data repository accession (to be added upon publication).

## Citation

If you use or reference this code, please cite the associated thesis/manuscript:

> Drelichman, M. *[Manuscript title — to be added upon publication].* San Diego State University.

## Author

**Maggie Drelichman** — San Diego State University
📧 mdrelichman@sdsu.edu

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
