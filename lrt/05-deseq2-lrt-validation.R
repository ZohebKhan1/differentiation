# 05-deseq2-lrt-validation.R
# Diagnostic visualizations for CPM+fold-filtered LRT gene lists from 03.
# Uses merged LS1+DSHC1 Di21 data for both heatmaps and CPM lineplots.
#
# Heatmaps: batch-corrected VST (limma), z-scored per gene.
# CPM lineplots: batch-corrected TMM CPM, x-axis split by experiment
#   at shared timepoints (D5, D7).
#
# Inputs:
#   results/LS1-DSHC1-merged/LRT_cpm_fold_filtered/{subset}/...csv
#   data/LS1-DSHC1-merged/merged_metadata.rds
#   data/LS1-DSHC1-merged/merged_raw_counts.rds
#
# Outputs (all in figures/LS1/maturation-score/validation/):
#   heatmap_cpmfold_{subset}.pdf       (x8)
#   top20_cpmfold_{subset}.svg         (x8)
#   bottom20_cpmfold_{subset}.svg      (x8, if enough genes)

library(DESeq2)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(grid)

# -- Parameters ----------------------------------------------------------------

page <- 1
genes_per_page <- 20

filtered_dir <- "results/LS1-DSHC1-merged/LRT_cpm_fold_filtered"
merged_data_dir <- "data/LS1-DSHC1-merged"
out_dir <- "figures/LS1/maturation-score/validation"

lrt_subsets <- c(
  "D3_D5_D7", "D3_D5_D7_D9", "D5_D7", "D5_D7_D9",
  "D5_D7_D9_D11", "D5_D7_D9_D11_D13", "D7_D9_D11", "D7_D9_D11_D13"
)

padj_map <- c(
  "D3_D5_D7" = 0.01, "D3_D5_D7_D9" = 0.01,
  "D5_D7" = 0.01, "D5_D7_D9" = 0.01,
  "D5_D7_D9_D11" = 0.001, "D5_D7_D9_D11_D13" = 0.001,
  "D7_D9_D11" = 0.001, "D7_D9_D11_D13" = 0.001
)

tp_colors <- c(
  "D3" = "#984EA3", "D5" = "#4DAF4A",
  "D7" = "red3", "D9" = "#FF7F00",
  "D11" = "#A65628", "D13" = "#E7298A"
)
experiment_colors <- c("LS1" = "#66C2A5", "DSHC1" = "#FC8D62")

# -- Load merged data ----------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

merged_meta <- readRDS(file.path(merged_data_dir, "merged_metadata.rds"))
merged_counts <- readRDS(file.path(merged_data_dir, "merged_raw_counts.rds"))

merged_meta$timepoint <- factor(
  merged_meta$timepoint,
  levels = c("D3", "D5", "D7", "D9", "D11", "D13")
)
merged_meta$experiment <- factor(
  merged_meta$experiment,
  levels = c("LS1", "DSHC1")
)

col_order <- order(merged_meta$timepoint, merged_meta$experiment)
merged_meta <- merged_meta[col_order, ]
merged_counts <- merged_counts[, merged_meta$sampleid, drop = FALSE]

stopifnot(all(colnames(merged_counts) == merged_meta$sampleid))

# -- Batch-corrected CPM (for lineplots) ---------------------------------------

dge <- edgeR::DGEList(counts = merged_counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")

log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
design_protect <- stats::model.matrix(~timepoint, data = merged_meta)
log_cpm <- limma::removeBatchEffect(
  log_cpm,
  batch = merged_meta$experiment,
  design = design_protect
)
merged_cpm <- 2^log_cpm - 1
merged_cpm[merged_cpm < 0] <- 0

# -- Batch-corrected VST (for heatmaps) ---------------------------------------

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = merged_counts,
  colData = merged_meta,
  design = ~ experiment + timepoint
)
vst_obj <- DESeq2::vst(dds, blind = FALSE)
vst_mat <- SummarizedExperiment::assay(vst_obj)

design_protect <- stats::model.matrix(~timepoint, data = merged_meta)
vst_corrected <- limma::removeBatchEffect(
  vst_mat,
  batch = merged_meta$experiment,
  design = design_protect
)

# -- Lineplot x-axis groups ----------------------------------------------------
# D3 = LS1 only; D9/D11/D13 = DSHC1 only; D5/D7 = both experiments

x_group <- ifelse(
  merged_meta$timepoint %in% c("D3", "D9", "D11", "D13"),
  as.character(merged_meta$timepoint),
  paste0(merged_meta$timepoint, "-", merged_meta$experiment)
)

x_levels <- c(
  "D3", "D5-LS1", "D5-DSHC1", "D7-LS1", "D7-DSHC1",
  "D9", "D11", "D13"
)
x_group <- factor(x_group, levels = x_levels)

# -- Heatmap setup -------------------------------------------------------------

heatmap_palette <- grDevices::colorRampPalette(c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
))(100)
heatmap_col <- circlize::colorRamp2(
  seq(-2, 2, length.out = 100), heatmap_palette
)

gpar_title <- grid::gpar(
  fontsize = 10, fontface = "bold", fontfamily = "sans"
)
gpar_labels <- grid::gpar(fontsize = 9, fontfamily = "sans")

col_anno <- ComplexHeatmap::HeatmapAnnotation(
  Timepoint = merged_meta$timepoint,
  Experiment = merged_meta$experiment,
  col = list(Timepoint = tp_colors, Experiment = experiment_colors),
  show_annotation_name = TRUE,
  annotation_name_gp = grid::gpar(
    fontsize = 9, fontfamily = "sans", col = "black"
  ),
  annotation_legend_param = list(
    Timepoint = list(
      title_gp = gpar_title, labels_gp = gpar_labels
    ),
    Experiment = list(
      title_gp = gpar_title, labels_gp = gpar_labels
    )
  )
)

col_split <- merged_meta$timepoint

# -- Lineplot builder ----------------------------------------------------------

lineplot_theme <- ggplot2::theme_bw(base_family = "sans", base_size = 8) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(
      size = 7, face = "bold", color = "black"
    ),
    strip.background = ggplot2::element_rect(fill = "grey95"),
    axis.title = ggplot2::element_text(size = 9, color = "black"),
    axis.text = ggplot2::element_text(size = 7, color = "black"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      size = 10, face = "bold", color = "black", hjust = 0
    ),
    plot.margin = ggplot2::margin(5, 5, 5, 5),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 8, color = "black"),
    legend.text = ggplot2::element_text(size = 7, color = "black")
  )

make_lineplot_grid <- function(gene_ids, plot_title) {
  plot_data <- do.call(rbind, lapply(gene_ids, function(g) {
    data.frame(
      gene = g,
      x_group = x_group,
      experiment = as.character(merged_meta$experiment),
      cpm = as.numeric(merged_cpm[g, ]),
      stringsAsFactors = FALSE
    )
  }))
  plot_data$gene <- factor(plot_data$gene, levels = gene_ids)
  plot_data$x_group <- factor(plot_data$x_group, levels = x_levels)
  plot_data$x_pos <- as.integer(plot_data$x_group)

  mean_data <- stats::aggregate(
    cpm ~ gene + x_group + experiment,
    data = plot_data, FUN = mean
  )
  mean_data$gene <- factor(mean_data$gene, levels = gene_ids)
  mean_data$x_group <- factor(mean_data$x_group, levels = x_levels)
  mean_data$x_pos <- as.integer(mean_data$x_group)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = x_pos, y = cpm)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = experiment),
      size = 1.2, alpha = 0.7
    ) +
    ggplot2::geom_line(
      data = mean_data,
      ggplot2::aes(
        x = x_pos, y = cpm,
        group = experiment, color = experiment
      ),
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = mean_data,
      ggplot2::aes(x = x_pos, y = cpm, color = experiment),
      size = 1.8
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(x_levels), labels = x_levels
    ) +
    ggplot2::scale_color_manual(
      values = experiment_colors, name = "Experiment"
    ) +
    ggplot2::facet_wrap(
      ~gene,
      nrow = 4, ncol = 5, scales = "free_y"
    ) +
    ggplot2::labs(x = "Timepoint", y = "CPM", title = plot_title) +
    lineplot_theme
}

# -- Generate heatmaps + lineplots per subset ----------------------------------

for (subset_name in lrt_subsets) {
  padj_threshold <- padj_map[[subset_name]]
  padj_label <- sub("^0\\.", "0.", format(padj_threshold, scientific = FALSE))

  gene_file <- file.path(
    filtered_dir, subset_name,
    paste0("padj", padj_label, "_genes_cpmFoldFiltered.csv")
  )
  sig_df <- read.csv(gene_file, stringsAsFactors = FALSE)

  # Heatmap
  genes_in_vst <- intersect(sig_df$gene_id, rownames(vst_corrected))
  mat_z <- t(apply(
    vst_corrected[genes_in_vst, , drop = FALSE], 1,
    function(row) {
      (row - mean(row, na.rm = TRUE)) / sd(row, na.rm = TRUE)
    }
  ))
  mat_z[mat_z < -2] <- -2
  mat_z[mat_z > 2] <- 2

  n_genes <- nrow(mat_z)
  pdf_height <- min(max(n_genes * 0.003 + 2, 6), 8)

  ht <- ComplexHeatmap::Heatmap(
    mat_z,
    name = "Z-score",
    col = heatmap_col,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    clustering_method_rows = "ward.D2",
    row_dend_reorder = TRUE,
    column_split = col_split,
    column_title = NULL,
    row_title = NULL,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = col_anno,
    show_row_dend = FALSE,
    border = TRUE,
    column_gap = grid::unit(1, "mm"),
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar_title,
      labels_gp = gpar_labels,
      at = c(-2, -1, 0, 1, 2),
      legend_height = grid::unit(3, "cm"),
      legend_direction = "vertical"
    ),
    use_raster = TRUE,
    raster_quality = 5
  )

  pdf(
    file.path(out_dir, paste0("heatmap_cpmfold_", subset_name, ".pdf")),
    width = 7, height = pdf_height
  )
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    padding = grid::unit(c(2, 8, 2, 2), "mm")
  )
  grDevices::dev.off()

  # CPM lineplots
  sig_genes <- sig_df[sig_df$gene_id %in% rownames(merged_cpm), ]
  n_sig <- nrow(sig_genes)
  if (n_sig == 0) next

  top_start <- (page - 1) * genes_per_page + 1
  top_end <- min(page * genes_per_page, n_sig)
  bottom_end <- n_sig - (page - 1) * genes_per_page
  bottom_start <- max(bottom_end - genes_per_page + 1, 1)

  if (top_start <= n_sig) {
    top_genes <- sig_genes$gene_id[top_start:top_end]
    p_top <- make_lineplot_grid(
      top_genes,
      sprintf(
        "%s cpmFold: Top genes %d-%d (padj<%s)",
        subset_name, top_start, top_end, padj_label
      )
    )
    ggplot2::ggsave(
      file.path(out_dir, paste0("top20_cpmfold_", subset_name, ".svg")),
      plot = p_top, device = svglite::svglite,
      width = 15, height = 8, fix_text_size = FALSE
    )
  }

  if (bottom_start >= 1 && bottom_start != top_start) {
    bottom_genes <- sig_genes$gene_id[bottom_start:bottom_end]
    p_bottom <- make_lineplot_grid(
      bottom_genes,
      sprintf(
        "%s cpmFold: Bottom genes %d-%d (padj<%s)",
        subset_name, bottom_start, bottom_end, padj_label
      )
    )
    ggplot2::ggsave(
      file.path(out_dir, paste0("bottom20_cpmfold_", subset_name, ".svg")),
      plot = p_bottom, device = svglite::svglite,
      width = 15, height = 8, fix_text_size = FALSE
    )
  }
}
