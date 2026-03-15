# 02-deseq2-lrt-LS1-DSHC1.R
# DESeq2 Likelihood Ratio Test on each CPM-filtered subset from 01.
# Tests whether timepoint explains variance beyond experiment batch.
#
# Design:  full = ~ experiment + timepoint, reduced = ~ experiment
#          (falls back to ~ timepoint / ~ 1 if only one experiment present)
#
# Inputs (from 01-filter-merge-LS1-DSHC1.R):
#   data/LS1-DSHC1-merged/subsets/{subset}/metadata.rds
#   data/LS1-DSHC1-merged/subsets/{subset}/raw_counts.rds
#
# Outputs (all in results/LS1-DSHC1-merged/LRT/):
#   {subset}/lrt_full_results.csv              -- all genes with LRT statistics
#   {subset}/padj0.05_genes.csv                -- genes with padj < 0.05
#   {subset}/padj0.01_genes.csv                -- genes with padj < 0.01
#   {subset}/padj0.001_genes.csv               -- genes with padj < 0.001
#   lrt_summary.csv                            -- gene counts across all subsets

source("scripts/LS1/preprocessing/00-dependencies-LS1.R")

# -- Parameters ----------------------------------------------------------------

input_dir <- "data/LS1-DSHC1-merged/subsets"
lrt_dir <- "results/LS1-DSHC1-merged/LRT"

lrt_subsets <- c(
  "D3_D5_D7", "D3_D5_D7_D9", "D5_D7", "D5_D7_D9",
  "D5_D7_D9_D11", "D5_D7_D9_D11_D13", "D7_D9_D11", "D7_D9_D11_D13"
)

padj_thresholds <- c(0.05, 0.01, 0.001)

# -- LRT runner ----------------------------------------------------------------

dir.create(lrt_dir, recursive = TRUE, showWarnings = FALSE)

run_lrt <- function(subset_name) {
  meta <- readRDS(file.path(input_dir, subset_name, "metadata.rds"))
  counts <- readRDS(file.path(input_dir, subset_name, "raw_counts.rds"))

  counts <- counts[, meta$sampleid, drop = FALSE]
  stopifnot(all(colnames(counts) == meta$sampleid))

  meta$timepoint <- droplevels(factor(meta$timepoint))
  meta$experiment <- droplevels(factor(meta$experiment))

  multi_experiment <- length(unique(meta$experiment)) > 1
  if (multi_experiment) {
    design_full <- ~ experiment + timepoint
    design_reduced <- ~experiment
  } else {
    design_full <- ~timepoint
    design_reduced <- ~1
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts, colData = meta, design = design_full
  )
  dds <- DESeq2::DESeq(dds, test = "LRT", reduced = design_reduced)

  res_df <- as.data.frame(DESeq2::results(dds))
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
  res_df <- res_df[, c(
    "gene_id", "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj"
  )]

  subset_out <- file.path(lrt_dir, subset_name)
  dir.create(subset_out, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(res_df, file.path(subset_out, "lrt_full_results.csv"),
    row.names = FALSE
  )

  for (threshold in padj_thresholds) {
    sig <- res_df[!is.na(res_df$padj) & res_df$padj < threshold, ]
    label <- sub("^0\\.", "0.", format(threshold, scientific = FALSE))
    utils::write.csv(
      sig,
      file.path(subset_out, paste0("padj", label, "_genes.csv")),
      row.names = FALSE
    )
  }

  data.frame(
    subset = subset_name,
    n_samples = ncol(counts),
    n_genes_tested = nrow(res_df),
    n_with_padj = sum(!is.na(res_df$padj)),
    sig_padj_0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05),
    sig_padj_0.01 = sum(!is.na(res_df$padj) & res_df$padj < 0.01),
    sig_padj_0.001 = sum(!is.na(res_df$padj) & res_df$padj < 0.001),
    batch_in_model = multi_experiment,
    stringsAsFactors = FALSE
  )
}

# -- Run across all subsets ----------------------------------------------------

summary_df <- do.call(rbind, lapply(lrt_subsets, run_lrt))
rownames(summary_df) <- NULL

utils::write.csv(summary_df, file.path(lrt_dir, "lrt_summary.csv"),
  row.names = FALSE
)
