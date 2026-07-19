## deseq2_pipeline/scripts/deseq2_pipeline.R
##
## Generic two-stage DESeq2 pipeline:
##   Part I  - DESeq2 on raw filtered counts
##   Part II - DESeq2 on ComBat-seq batch-corrected counts (computed here, not
##             read from an external file, so the whole pipeline is reproducible
##             from raw counts alone)
##
## Edit the CONFIG block below to point at a new dataset -- nothing else in
## this script should need to change. Paths are relative to the repo root
## (i.e. run this with deseq2_pipeline/ as your working directory's parent).

## ---- CONFIG ---------------------------------------------------------------

counts_path   <- "deseq2_pipeline/data/2023.GFP18_AKT_FeatureCounts.csv"
coldata_path  <- "deseq2_pipeline/data/2023.GFP18_AKT_annotations.csv"

condition_col <- "Condition"   # colData column DESeq2's design is built on
reference_level <- "control"  # baseline level for condition_col
batch_col     <- "Batch"      # colData column ComBat-seq corrects on -- set to
                               # whatever batch variable exists in your colData

min_count_sum <- 10            # low-count filter: keep genes with rowSums(counts) > this
alpha         <- 0.05           # significance threshold for results()

out_dir <- "deseq2_pipeline/results"   # all tables/plots are written under here

## ---- Packages ---------------------------------------------------------

required_bioc <- c("DESeq2", "sva")
required_cran <- c("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
}
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
  library(ggplot2)
})

## ---- Load + align data -----------------------------------------------

counts_data <- read.csv(counts_path, header = TRUE, row.names = 1, check.names = FALSE)
colData     <- read.csv(coldata_path, header = TRUE, row.names = 1)

counts_data <- counts_data[rowSums(counts_data) > min_count_sum, ]
message(sprintf("Kept %d genes after low-count filter (rowSums > %d)",
                nrow(counts_data), min_count_sum))

## DESeqDataSetFromMatrix does NOT align samples for you -- it assumes
## colData row i corresponds to counts column i. Enforce that explicitly
## instead of just checking set membership.
stopifnot("Every count-matrix sample must have a colData row" =
            all(colnames(counts_data) %in% rownames(colData)))
colData <- colData[colnames(counts_data), , drop = FALSE]
stopifnot("colData rows and counts columns are out of order" =
            all(rownames(colData) == colnames(counts_data)))

colData[[condition_col]] <- relevel(factor(colData[[condition_col]]), ref = reference_level)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Reusable DESeq2 runner --------------------------------------------

run_deseq2 <- function(counts, colData, design, alpha, label) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = design)
  dds <- DESeq(dds)
  res <- results(dds, alpha = alpha)

  message(sprintf("--- %s: DESeq2 summary (alpha = %.2f) ---", label, alpha))
  print(summary(res))

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  write.csv(res_df, file.path(out_dir, sprintf("%s_results.csv", label)), row.names = FALSE)

  png(file.path(out_dir, sprintf("%s_MA_plot.png", label)), width = 800, height = 600)
  plotMA(res, main = sprintf("%s: MA plot", label))
  dev.off()

  vsd <- vst(dds, blind = TRUE)
  pca <- plotPCA(vsd, intgroup = condition_col)
  ggsave(file.path(out_dir, sprintf("%s_PCA.png", label)), pca, width = 6, height = 5, dpi = 150)

  list(dds = dds, res = res)
}

## ---- Part I: DESeq2 on raw (uncorrected) counts ------------------------

design_formula <- as.formula(paste("~", condition_col))
part1 <- run_deseq2(counts_data, colData, design_formula, alpha, label = "raw")

## ---- Part II: ComBat-seq batch correction, then DESeq2 -----------------

stopifnot("batch_col must exist in colData to run ComBat-seq" = batch_col %in% colnames(colData))

combat_counts <- ComBat_seq(
  counts = as.matrix(counts_data),
  batch  = colData[[batch_col]],
  group  = colData[[condition_col]]
)
write.csv(combat_counts, file.path(out_dir, "combat_seq_corrected_counts.csv"))

part2 <- run_deseq2(combat_counts, colData, design_formula, alpha, label = "combat_corrected")

## ---- Session info (for reproducibility) --------------------------------

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message(sprintf("Done. Results, plots, and sessionInfo written to %s/", out_dir))
