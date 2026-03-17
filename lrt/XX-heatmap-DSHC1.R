# XX-heatmap-DSHC1.R
# heatmap of LS1 Di21 maturation genes (D7 vs D5,
# log2FC > 0.5, padj < 0.05) plotted in DSHC1
# Di21 samples (D5, D7, D9).
#
# VST row z-scored, winsorized [-2, 2],
# quantile-mapped blue-red diverging.
# columns: D5 ctrl | D7 ctrl + SAG | D9 ctrl + SAG
#
# inputs:
#   DSHC1_vst, DSHC1_metadata (from loader)
#   data/LS1/degs/di21-late-vs-early/
#     LS1_degs_Di21ctrlD7vDi21ctrlD5_filtered.csv
#
# outputs:
#   figures/DSHC1/heatmap/
#     LS1_maturation_D5toD7_in_DSHC1_heatmap.png
#     LS1_progenitor_D5toD7_in_DSHC1_heatmap.png

source(
  'scripts/DSHC1/00-preprocessing/00-load-data-DSHC1.R'
)
sysfonts::font_add(
  'Arial',
  regular = 'C:/Windows/Fonts/arial.ttf',
  bold = 'C:/Windows/Fonts/arialbd.ttf',
  italic = 'C:/Windows/Fonts/ariali.ttf'
)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)

# --- parameters ----

z_cap <- 2

# figure dimensions
fig_width <- 2.953
fig_height <- 3.74

# fonts
font_family <- 'Arial'
tp_label_size <- 8
legend_title_size <- 8
legend_text_size <- 8

# column gaps (mm)
between_tp_gap <- 3
within_tp_gap <- 0.75

# output path
fig_dir <- 'figures/DSHC1/heatmap'
dir.create(
  fig_dir,
  showWarnings = FALSE, recursive = TRUE
)

# --- load LS1 maturation genes ----

ls1_deg_path <- paste0(
  'data/LS1/degs/di21-late-vs-early/',
  'LS1_degs_Di21ctrlD7vDi21ctrlD5_filtered.csv'
)
ls1_degs <- read.csv(ls1_deg_path)

# maturation = upregulated D5 -> D7
maturation_genes <- ls1_degs$gene_id[
  ls1_degs$log2FoldChange > 0.5 &
    ls1_degs$padj < 0.05
]

# intersect with DSHC1 VST rownames
maturation_genes <- intersect(
  maturation_genes, rownames(DSHC1_vst)
)

# progenitor = downregulated D5 -> D7
progenitor_genes <- ls1_degs$gene_id[
  ls1_degs$log2FoldChange < -0.5 &
    ls1_degs$padj < 0.05
]
progenitor_genes <- intersect(
  progenitor_genes, rownames(DSHC1_vst)
)

# --- subset DSHC1 metadata for Di21 ----

di21_meta <- DSHC1_metadata[
  DSHC1_metadata$genotype == 'Di21', ,
  drop = FALSE
]

ctrl_d5 <- di21_meta[
  di21_meta$treatment == 'ctrl' &
    di21_meta$timepoint == 'D5', ,
  drop = FALSE
]
ctrl_d7 <- di21_meta[
  di21_meta$treatment == 'ctrl' &
    di21_meta$timepoint == 'D7', ,
  drop = FALSE
]
sag_d7 <- di21_meta[
  di21_meta$treatment == 'SAG' &
    di21_meta$timepoint == 'D7', ,
  drop = FALSE
]
ctrl_d9 <- di21_meta[
  di21_meta$treatment == 'ctrl' &
    di21_meta$timepoint == 'D9', ,
  drop = FALSE
]
sag_d9 <- di21_meta[
  di21_meta$treatment == 'SAG' &
    di21_meta$timepoint == 'D9', ,
  drop = FALSE
]

plot_meta <- rbind(
  ctrl_d5, ctrl_d7, sag_d7, ctrl_d9, sag_d9
)
plot_samples <- rownames(plot_meta)

# --- build z-scored matrix ----

vst_sub <- DSHC1_vst[maturation_genes, plot_samples]

vst_z <- t(apply(vst_sub, 1, function(row) {
  mu <- mean(row, na.rm = TRUE)
  s <- sd(row, na.rm = TRUE)
  if (s == 0) {
    return(rep(0, length(row)))
  }
  (row - mu) / s
}))

# winsorize
vst_z[vst_z < -z_cap] <- -z_cap
vst_z[vst_z > z_cap] <- z_cap

# --- column split ----

col_split <- factor(
  paste0(
    plot_meta$treatment, '_',
    plot_meta$timepoint
  ),
  levels = c(
    'ctrl_D5', 'ctrl_D7', 'SAG_D7',
    'ctrl_D9', 'SAG_D9'
  )
)

# 4 gaps: big, small, big, small
col_gaps <- grid::unit(
  c(
    between_tp_gap, within_tp_gap,
    between_tp_gap, within_tp_gap
  ),
  'mm'
)

# --- color scale (quantile blue-red) ----

all_values <- as.vector(vst_z)

quantile_breaks <- quantile(
  all_values,
  probs = seq(0, 1, length.out = 100),
  na.rm = TRUE
)

blue_red_palette <- grDevices::colorRampPalette(c(
  '#053061', '#2166AC', '#4393C3',
  '#92C5DE', '#D1E5F0',
  '#FDDBC7', '#F4A582', '#D6604D',
  '#B2182B', '#67001F'
))(100)

col_fun <- circlize::colorRamp2(
  quantile_breaks, blue_red_palette
)

# --- build heatmap ----

ht <- ComplexHeatmap::Heatmap(
  vst_z,
  name = 'Z-score',
  col = col_fun,
  column_split = col_split,
  cluster_rows = TRUE,
  clustering_method_rows = 'ward.D2',
  cluster_columns = FALSE,
  column_title = c(
    'D5\nctrl', 'D7\nctrl', 'D7\nSAG',
    'D9\nctrl', 'D9\nSAG'
  ),
  column_title_side = 'bottom',
  column_title_gp = grid::gpar(
    fontsize = tp_label_size,
    fontfamily = font_family,
    col = 'black'
  ),
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  row_dend_side = 'left',
  row_dend_reorder = TRUE,
  row_dend_width = grid::unit(8, 'mm'),
  border = FALSE,
  rect_gp = grid::gpar(col = NA),
  column_gap = col_gaps,
  heatmap_legend_param = list(
    title = 'Z-score',
    title_gp = grid::gpar(
      fontsize = legend_title_size,
      fontfamily = font_family,
      col = 'black'
    ),
    labels_gp = grid::gpar(
      fontsize = legend_text_size,
      fontfamily = font_family,
      col = 'black'
    ),
    at = c(-z_cap, -1, 0, 1, z_cap),
    labels = c(
      paste0('\u2264', '-', z_cap),
      '-1', '0', '1',
      paste0('\u2265', z_cap)
    ),
    legend_height = grid::unit(3, 'cm'),
    legend_width = grid::unit(3, 'mm'),
    legend_direction = 'vertical'
  ),
  use_raster = TRUE,
  raster_quality = 5
)

# --- save as high-res PNG ----

output_fname <- 'LS1_maturation_D5toD7_in_DSHC1_heatmap.png'

grDevices::png(
  file.path(fig_dir, output_fname),
  width = fig_width,
  height = fig_height,
  units = 'in',
  res = 600
)
ComplexHeatmap::draw(
  ht,
  show_annotation_legend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_side = 'right',
  padding = grid::unit(c(2, 10, 2, 2), 'mm')
)
grDevices::dev.off()

# =============================================================
# progenitor heatmap (decreasing D5 -> D7)
# =============================================================

vst_sub_prog <- DSHC1_vst[
  progenitor_genes, plot_samples
]

vst_z_prog <- t(apply(
  vst_sub_prog, 1, function(row) {
    mu <- mean(row, na.rm = TRUE)
    s <- sd(row, na.rm = TRUE)
    if (s == 0) {
      return(rep(0, length(row)))
    }
    (row - mu) / s
  }
))

vst_z_prog[vst_z_prog < -z_cap] <- -z_cap
vst_z_prog[vst_z_prog > z_cap] <- z_cap

# quantile color scale for progenitor
all_vals_prog <- as.vector(vst_z_prog)
qbreaks_prog <- quantile(
  all_vals_prog,
  probs = seq(0, 1, length.out = 100),
  na.rm = TRUE
)
col_fun_prog <- circlize::colorRamp2(
  qbreaks_prog, blue_red_palette
)

ht_prog <- ComplexHeatmap::Heatmap(
  vst_z_prog,
  name = 'Z-score',
  col = col_fun_prog,
  column_split = col_split,
  cluster_rows = TRUE,
  clustering_method_rows = 'ward.D2',
  cluster_columns = FALSE,
  column_title = c(
    'D5\nctrl', 'D7\nctrl', 'D7\nSAG',
    'D9\nctrl', 'D9\nSAG'
  ),
  column_title_side = 'bottom',
  column_title_gp = grid::gpar(
    fontsize = tp_label_size,
    fontfamily = font_family,
    col = 'black'
  ),
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  row_dend_side = 'left',
  row_dend_reorder = TRUE,
  row_dend_width = grid::unit(8, 'mm'),
  border = FALSE,
  rect_gp = grid::gpar(col = NA),
  column_gap = col_gaps,
  heatmap_legend_param = list(
    title = 'Z-score',
    title_gp = grid::gpar(
      fontsize = legend_title_size,
      fontfamily = font_family,
      col = 'black'
    ),
    labels_gp = grid::gpar(
      fontsize = legend_text_size,
      fontfamily = font_family,
      col = 'black'
    ),
    at = c(-z_cap, -1, 0, 1, z_cap),
    labels = c(
      paste0('\u2264', '-', z_cap),
      '-1', '0', '1',
      paste0('\u2265', z_cap)
    ),
    legend_height = grid::unit(3, 'cm'),
    legend_width = grid::unit(3, 'mm'),
    legend_direction = 'vertical'
  ),
  use_raster = TRUE,
  raster_quality = 5
)

prog_fname <- 'LS1_progenitor_D5toD7_in_DSHC1_heatmap.png'

grDevices::png(
  file.path(fig_dir, prog_fname),
  width = fig_width,
  height = fig_height,
  units = 'in',
  res = 600
)
ComplexHeatmap::draw(
  ht_prog,
  show_annotation_legend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_side = 'right',
  padding = grid::unit(c(2, 10, 2, 2), 'mm')
)
grDevices::dev.off()

# =============================================================
# SAG-dependent DEG heatmaps (combined induced + repressed)
# Di21 SAG vs Di21 ctrl, row_split by direction
# =============================================================

sag_deg_dir <- 'data/DSHC1/degs/di21-sag-vs-di21-ctrl'
sag_row_gap_mm <- 3

make_sag_combined_heatmap <- function(
  deg_csv_path, output_name
) {
  degs <- read.csv(deg_csv_path)

  up_genes <- degs$gene_id[
    degs$log2FoldChange > 0.5 &
      degs$padj < 0.05
  ]
  down_genes <- degs$gene_id[
    degs$log2FoldChange < -0.5 &
      degs$padj < 0.05
  ]

  up_genes <- intersect(
    up_genes, rownames(DSHC1_vst)
  )
  down_genes <- intersect(
    down_genes, rownames(DSHC1_vst)
  )

  all_genes <- c(up_genes, down_genes)

  mat <- DSHC1_vst[all_genes, plot_samples]
  mat_z <- t(apply(mat, 1, function(row) {
    mu <- mean(row, na.rm = TRUE)
    s <- sd(row, na.rm = TRUE)
    if (s == 0) {
      return(rep(0, length(row)))
    }
    (row - mu) / s
  }))
  mat_z[mat_z < -z_cap] <- -z_cap
  mat_z[mat_z > z_cap] <- z_cap

  row_split <- factor(
    c(
      rep('SAG-induced', length(up_genes)),
      rep('SAG-repressed', length(down_genes))
    ),
    levels = c('SAG-induced', 'SAG-repressed')
  )

  vals <- as.vector(mat_z)
  qb <- quantile(
    vals,
    probs = seq(0, 1, length.out = 100),
    na.rm = TRUE
  )
  cf <- circlize::colorRamp2(qb, blue_red_palette)

  ht_sag <- ComplexHeatmap::Heatmap(
    mat_z,
    name = 'Z-score',
    col = cf,
    row_split = row_split,
    column_split = col_split,
    cluster_rows = TRUE,
    clustering_method_rows = 'ward.D2',
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    column_title = c(
      'D5\nctrl', 'D7\nctrl', 'D7\nSAG',
      'D9\nctrl', 'D9\nSAG'
    ),
    column_title_side = 'bottom',
    column_title_gp = grid::gpar(
      fontsize = tp_label_size,
      fontfamily = font_family,
      col = 'black'
    ),
    row_title = c('SAG-induced', 'SAG-repressed'),
    row_title_gp = grid::gpar(
      fontsize = 8,
      fontfamily = font_family,
      col = 'black'
    ),
    row_title_rot = 270,
    row_title_side = 'right',
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_heatmap_legend = TRUE,
    row_dend_side = 'left',
    row_dend_reorder = TRUE,
    row_dend_width = grid::unit(8, 'mm'),
    border = FALSE,
    rect_gp = grid::gpar(col = NA),
    column_gap = col_gaps,
    row_gap = grid::unit(sag_row_gap_mm, 'mm'),
    heatmap_legend_param = list(
      title = 'Z-score',
      title_gp = grid::gpar(
        fontsize = legend_title_size,
        fontfamily = font_family,
        col = 'black'
      ),
      labels_gp = grid::gpar(
        fontsize = legend_text_size,
        fontfamily = font_family,
        col = 'black'
      ),
      at = c(-z_cap, -1, 0, 1, z_cap),
      labels = c(
        paste0('\u2264', '-', z_cap),
        '-1', '0', '1',
        paste0('\u2265', z_cap)
      ),
      legend_height = grid::unit(3, 'cm'),
      legend_width = grid::unit(3, 'mm'),
      legend_direction = 'vertical'
    ),
    use_raster = TRUE,
    raster_quality = 5
  )

  grDevices::png(
    file.path(fig_dir, output_name),
    width = fig_width,
    height = fig_height,
    units = 'in',
    res = 600
  )
  ComplexHeatmap::draw(
    ht_sag,
    show_annotation_legend = FALSE,
    show_heatmap_legend = TRUE,
    heatmap_legend_side = 'right',
    padding = grid::unit(c(2, 10, 2, 2), 'mm')
  )
  grDevices::dev.off()
}

# D7: SAG induced + repressed combined
make_sag_combined_heatmap(
  file.path(
    sag_deg_dir,
    'Di21_SAG_D7_vs_Di21_ctrl_D7_filtered.csv'
  ),
  'Di21_SAG_D7_heatmap.png'
)

# D9: SAG induced + repressed combined
make_sag_combined_heatmap(
  file.path(
    sag_deg_dir,
    'Di21_SAG_D9_vs_Di21_ctrl_D9_filtered.csv'
  ),
  'Di21_SAG_D9_heatmap.png'
)
