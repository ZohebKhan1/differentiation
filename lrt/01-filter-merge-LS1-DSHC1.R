# 01-filter-merge-LS1-DSHC1.R
# Merge LS1 + DSHC1 Di21 control samples and create CPM-filtered subsets
# for downstream LRT analysis.
#
# Inputs:
#   scripts/LS1/preprocessing/00-load-data-LS1.R -- LS1 metadata + raw counts
#   data/DSHC1/metadata/DSHC1_metadata.rds      -- DSHC1 sample metadata
#   data/DSHC1/DSHC1_dds.rds                    -- DSHC1 fitted DESeqDataSet
#
# Outputs (all in data/LS1-DSHC1-merged/):
#   merged_metadata.rds                         -- 24-sample unified metadata
#   merged_raw_counts.rds                       -- 24-sample raw count matrix
#   subsets/{subset_name}/metadata.rds           -- per-subset metadata (x8)
#   subsets/{subset_name}/raw_counts.rds         -- per-subset filtered counts (x8)

source("scripts/LS1/preprocessing/00-dependencies-LS1.R")
source("scripts/LS1/preprocessing/00-load-data-LS1.R")

# -- Parameters ----------------------------------------------------------------

out_dir <- "data/LS1-DSHC1-merged"
cpm_threshold <- 5

lrt_subsets <- list(
  "D3_D5_D7"         = c("D3", "D5", "D7"),
  "D3_D5_D7_D9"      = c("D3", "D5", "D7", "D9"),
  "D5_D7"            = c("D5", "D7"),
  "D5_D7_D9"         = c("D5", "D7", "D9"),
  "D5_D7_D9_D11"     = c("D5", "D7", "D9", "D11"),
  "D5_D7_D9_D11_D13" = c("D5", "D7", "D9", "D11", "D13"),
  "D7_D9_D11"        = c("D7", "D9", "D11"),
  "D7_D9_D11_D13"    = c("D7", "D9", "D11", "D13")
)

# -- Load LS1 (Experiment 1, Set 3) -- Di21 at D3/D5/D7 -----------------------

ls1_meta <- LS1_metadata
ls1_raw <- LS1_raw_counts

ls1_keep <- ls1_meta$genotype == "Di21" &
  ls1_meta$timepoint %in% c("D3", "D5", "D7")
ls1_meta <- ls1_meta[ls1_keep, ]
ls1_raw <- ls1_raw[, ls1_meta$sampleid, drop = FALSE]

stopifnot(nrow(ls1_meta) == 9L)
stopifnot(all(colnames(ls1_raw) == ls1_meta$sampleid))

# -- Load DSHC1 (Experiment 2) -- Di21 ctrl at D5/D7/D9/D11/D13 --------------

dshc1_meta <- readRDS("data/DSHC1/metadata/DSHC1_metadata.rds")
dshc1_raw <- DESeq2::counts(
  readRDS("data/DSHC1/DSHC1_dds.rds"),
  normalized = FALSE
)

dshc1_keep <- dshc1_meta$genotype == "Di21" &
  dshc1_meta$treatment == "ctrl" &
  dshc1_meta$timepoint %in% c("D5", "D7", "D9", "D11", "D13")
dshc1_meta <- dshc1_meta[dshc1_keep, ]
dshc1_raw <- dshc1_raw[, dshc1_meta$sampleid, drop = FALSE]

stopifnot(nrow(dshc1_meta) == 15L)
stopifnot(all(colnames(dshc1_raw) == dshc1_meta$sampleid))

# -- Harmonize metadata -------------------------------------------------------

make_unified <- function(meta) {
  data.frame(
    sampleid             = meta$sampleid,
    experiment           = meta$experiment_id,
    timepoint            = meta$timepoint,
    timepoint_numeric    = meta$timepoint_numeric,
    biological_replicate = meta$biological_replicate,
    stringsAsFactors     = FALSE
  )
}

merged_meta <- rbind(make_unified(ls1_meta), make_unified(dshc1_meta))
rownames(merged_meta) <- merged_meta$sampleid
merged_meta$experiment <- factor(merged_meta$experiment, levels = c("LS1", "DSHC1"))
merged_meta$timepoint <- factor(
  merged_meta$timepoint,
  levels = c("D3", "D5", "D7", "D9", "D11", "D13")
)

# -- Merge counts on common genes ---------------------------------------------

common_genes <- intersect(rownames(ls1_raw), rownames(dshc1_raw))
merged_counts <- cbind(
  ls1_raw[common_genes, , drop = FALSE],
  dshc1_raw[common_genes, , drop = FALSE]
)
merged_counts <- merged_counts[, merged_meta$sampleid, drop = FALSE]

stopifnot(ncol(merged_counts) == 24L)
stopifnot(all(colnames(merged_counts) == merged_meta$sampleid))
stopifnot(!anyDuplicated(merged_meta$sampleid))
# LS1 counts are prefiltered in preprocessing, so merged gene count is lower
# than legacy unfiltered matrices but should still be comfortably large.
stopifnot(nrow(merged_counts) > 5000L)

tp_counts <- table(merged_meta$timepoint)
stopifnot(tp_counts[["D3"]] == 3L)
stopifnot(tp_counts[["D5"]] == 6L)
stopifnot(tp_counts[["D7"]] == 6L)
stopifnot(tp_counts[["D9"]] == 3L)
stopifnot(tp_counts[["D11"]] == 3L)
stopifnot(tp_counts[["D13"]] == 3L)

# -- Per-subset CPM filtering -------------------------------------------------
# Each subset gets its own TMM normalization. A gene passes if
# meanCPM >= threshold at ANY timepoint within the subset.

dir.create(file.path(out_dir, "subsets"), recursive = TRUE, showWarnings = FALSE)

for (subset_name in names(lrt_subsets)) {
  subset_timepoints <- lrt_subsets[[subset_name]]

  subset_idx <- merged_meta$timepoint %in% subset_timepoints
  subset_meta <- merged_meta[subset_idx, ]
  subset_counts <- merged_counts[, subset_meta$sampleid, drop = FALSE]

  subset_meta$timepoint <- droplevels(subset_meta$timepoint)
  subset_meta$experiment <- droplevels(subset_meta$experiment)

  dge <- edgeR::DGEList(counts = subset_counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  cpm_mat <- edgeR::cpm(dge, log = FALSE)

  mean_cpm_by_tp <- sapply(levels(subset_meta$timepoint), function(tp) {
    tp_samples <- subset_meta$sampleid[subset_meta$timepoint == tp]
    rowMeans(cpm_mat[, tp_samples, drop = FALSE])
  })

  passes <- apply(mean_cpm_by_tp >= cpm_threshold, 1, any)
  filtered_counts <- subset_counts[passes, , drop = FALSE]

  if (nrow(filtered_counts) <= 5000L) {
    stop(sprintf(
      "%s: only %d genes passed CPM filter (expected > 5000)",
      subset_name, nrow(filtered_counts)
    ))
  }

  subset_out_dir <- file.path(out_dir, "subsets", subset_name)
  dir.create(subset_out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(subset_meta, file.path(subset_out_dir, "metadata.rds"))
  saveRDS(filtered_counts, file.path(subset_out_dir, "raw_counts.rds"))
}

# -- Save merged (unfiltered) data ---------------------------------------------

saveRDS(merged_meta, file.path(out_dir, "merged_metadata.rds"))
saveRDS(merged_counts, file.path(out_dir, "merged_raw_counts.rds"))
