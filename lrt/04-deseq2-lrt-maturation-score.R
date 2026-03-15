# 04-deseq2-lrt-maturation-score.R
# Score LS1 Di21 + Tri21 samples (D3/D5/D7) along maturation axes defined
# by CPM+fold-change filtered LRT gene lists from 03.
# PCA is trained on Di21 reference only; Tri21 is projected into that space.
# Two scoring methods: best-fit axis and centroid vector.
#
# Inputs:
#   results/LS1-DSHC1-merged/LRT_cpm_fold_filtered/{subset}/...csv
#   data/LS1/counts/LS1_vst_corrected.rds
#   data/LS1/metadata/LS1_metadata.rds
#
# Outputs:
#   results/LS1/maturation-score/maturation_scores.csv
#   figures/LS1/maturation-score/{subset}_{method}.svg      (x20)

library(ggplot2)

# -- Parameters ----------------------------------------------------------------

lrt_dir <- "results/LS1-DSHC1-merged/LRT_cpm_fold_filtered"
out_dir <- "results/LS1/maturation-score"
fig_dir <- "figures/LS1/maturation-score"

vst_path <- "data/LS1/counts/LS1_vst_corrected.rds"
meta_path <- "data/LS1/metadata/LS1_metadata.rds"

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

genotype_colors <- c("Di21" = "royalblue3", "Tri21" = "red3")

# Anchors for normalizing best-fit scores (0 = start, 1 = end)
bestfit_anchors <- list(
  "D3_D5_D7"         = c(start = 3L, end = 7L),
  "D3_D5_D7_D9"      = c(start = 3L, end = 7L),
  "D5_D7"            = c(start = 5L, end = 7L),
  "D5_D7_D9"         = c(start = 5L, end = 7L),
  "D5_D7_D9_D11"     = c(start = 5L, end = 7L),
  "D5_D7_D9_D11_D13" = c(start = 5L, end = 7L),
  "D7_D9_D11"        = c(start = 5L, end = 7L),
  "D7_D9_D11_D13"    = c(start = 5L, end = 7L)
)

# Centroid vector pairs (each list element = list of c(start_tp, end_tp))
centroid_vectors <- list(
  "D3_D5_D7"         = list(c(3L, 5L), c(5L, 7L), c(3L, 7L)),
  "D3_D5_D7_D9"      = list(c(3L, 5L), c(5L, 7L), c(3L, 7L)),
  "D5_D7"            = list(c(5L, 7L)),
  "D5_D7_D9"         = list(c(5L, 7L)),
  "D5_D7_D9_D11"     = list(c(5L, 7L)),
  "D5_D7_D9_D11_D13" = list(c(5L, 7L)),
  "D7_D9_D11"        = list(c(5L, 7L)),
  "D7_D9_D11_D13"    = list(c(5L, 7L))
)

# -- Helper functions ----------------------------------------------------------

flip_pc1 <- function(pca_fit, tp_numeric) {
  pc1 <- pca_fit$x[, "PC1"]
  m_early <- mean(pc1[tp_numeric == min(tp_numeric)])
  m_late <- mean(pc1[tp_numeric == max(tp_numeric)])
  if (m_late < m_early) {
    pca_fit$rotation[, "PC1"] <- -pca_fit$rotation[, "PC1"]
    pca_fit$x[, "PC1"] <- -pca_fit$x[, "PC1"]
  }
  pca_fit
}

project_samples <- function(vst_mat, pca_ref, n_pcs = 3L) {
  x_centered <- scale(t(vst_mat), center = pca_ref$center, scale = FALSE)
  coords <- x_centered %*% pca_ref$rotation[, seq_len(n_pcs), drop = FALSE]
  as.data.frame(coords)
}

bestfit_axis <- function(ref_coords) {
  x <- as.matrix(ref_coords[, c("PC1", "PC2", "PC3")])
  mu <- colMeans(x)
  xc <- sweep(x, 2, mu)
  v <- stats::prcomp(xc, center = FALSE, scale. = FALSE)$rotation[, 1]
  list(mu = mu, v = v)
}

centroid_axis <- function(ref_coords, start_tp, end_tp) {
  pcs <- c("PC1", "PC2", "PC3")
  mu1 <- colMeans(as.matrix(ref_coords[ref_coords$tp == start_tp, pcs]))
  mu2 <- colMeans(as.matrix(ref_coords[ref_coords$tp == end_tp, pcs]))
  v_raw <- mu2 - mu1
  l_dist <- sqrt(sum(v_raw^2))
  list(mu = mu1, v = v_raw / l_dist, l_dist = l_dist)
}

# -- Load data -----------------------------------------------------------------

vst_full <- readRDS(vst_path)
meta_full <- readRDS(meta_path)

meta <- meta_full[meta_full$timepoint %in% c("D3", "D5", "D7"), ]
meta <- droplevels(meta)
vst <- vst_full[, meta$sampleid, drop = FALSE]

ref_ids <- meta$sampleid[meta$genotype == "Di21"]
ref_meta <- meta[meta$genotype == "Di21", ]

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Scoring loop --------------------------------------------------------------

all_scores <- list()

for (lrt_sub in lrt_subsets) {
  padj_threshold <- padj_map[[lrt_sub]]
  padj_label <- sub("^0\\.", "0.", format(padj_threshold, scientific = FALSE))

  gene_file <- file.path(
    lrt_dir, lrt_sub,
    paste0("padj", padj_label, "_genes_cpmFoldFiltered.csv")
  )
  gene_df <- utils::read.csv(gene_file, stringsAsFactors = FALSE)
  genes <- intersect(gene_df$gene_id, rownames(vst))
  n_genes <- length(genes)

  pca_ref <- stats::prcomp(t(vst[genes, ref_ids, drop = FALSE]),
    center = TRUE, scale. = FALSE
  )
  pca_ref <- flip_pc1(pca_ref, ref_meta$timepoint_numeric)

  proj <- project_samples(vst[genes, , drop = FALSE], pca_ref, n_pcs = 3L)
  proj$tp <- meta$timepoint_numeric[match(rownames(proj), meta$sampleid)]

  ref_proj <- proj[ref_ids, ]

  # Best-fit scoring
  axis_bf <- bestfit_axis(ref_proj)
  t_raw <- as.numeric(sweep(
    as.matrix(proj[, c("PC1", "PC2", "PC3")]),
    2, axis_bf$mu
  ) %*% axis_bf$v)
  anchors <- bestfit_anchors[[lrt_sub]]
  is_ref <- rownames(proj) %in% ref_ids
  m0 <- mean(t_raw[is_ref & proj$tp == anchors[["start"]]])
  m1 <- mean(t_raw[is_ref & proj$tp == anchors[["end"]]])
  scores_bf <- (t_raw - m0) / (m1 - m0)

  all_scores[[length(all_scores) + 1L]] <- data.frame(
    sampleid = rownames(proj),
    genotype = meta$genotype[match(rownames(proj), meta$sampleid)],
    timepoint = meta$timepoint[match(rownames(proj), meta$sampleid)],
    timepoint_numeric = proj$tp,
    maturation_score = scores_bf,
    lrt_subset = lrt_sub,
    padj_threshold = padj_threshold,
    method = "best_fit",
    vector_label = "best_fit",
    n_genes = n_genes,
    stringsAsFactors = FALSE
  )

  # Centroid vector scoring
  for (vec in centroid_vectors[[lrt_sub]]) {
    axis_cv <- centroid_axis(ref_proj, start_tp = vec[1], end_tp = vec[2])
    t_raw_cv <- as.numeric(sweep(
      as.matrix(proj[, c("PC1", "PC2", "PC3")]),
      2, axis_cv$mu
    ) %*% axis_cv$v)
    scores_cv <- t_raw_cv / axis_cv$l_dist
    vec_label <- paste0("D", vec[1], "toD", vec[2])

    all_scores[[length(all_scores) + 1L]] <- data.frame(
      sampleid = rownames(proj),
      genotype = meta$genotype[match(rownames(proj), meta$sampleid)],
      timepoint = meta$timepoint[match(rownames(proj), meta$sampleid)],
      timepoint_numeric = proj$tp,
      maturation_score = scores_cv,
      lrt_subset = lrt_sub,
      padj_threshold = padj_threshold,
      method = "centroid",
      vector_label = vec_label,
      n_genes = n_genes,
      stringsAsFactors = FALSE
    )
  }
}

scores_df <- do.call(rbind, all_scores)
rownames(scores_df) <- NULL
utils::write.csv(scores_df, file.path(out_dir, "maturation_scores.csv"),
  row.names = FALSE
)

# -- Boxplot generation --------------------------------------------------------

make_panel <- function(df, title_str) {
  df$gene_label <- paste0("padj < ", df$padj_threshold, "  (n=", df$n_genes, ")")
  df$timepoint <- factor(df$timepoint, levels = c("D3", "D5", "D7"))
  df$genotype <- factor(df$genotype, levels = c("Di21", "Tri21"))

  ggplot2::ggplot(df, ggplot2::aes(
    x = genotype, y = maturation_score,
    fill = genotype
  )) +
    ggplot2::geom_boxplot(
      width = 0.5, outlier.shape = NA, linewidth = 0.4,
      colour = "black"
    ) +
    ggplot2::geom_point(
      shape = 21, size = 2, stroke = 0.3,
      colour = "black",
      position = ggplot2::position_jitter(width = 0.1, seed = 42)
    ) +
    ggplot2::facet_grid(timepoint ~ gene_label, scales = "free_y") +
    ggplot2::scale_fill_manual(values = genotype_colors) +
    ggplot2::labs(title = title_str, x = NULL, y = "Maturation Score") +
    ggplot2::theme_bw(base_family = "Arial", base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black", face = "bold", size = 11
      ),
      axis.title = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(color = "black"),
      legend.position = "none"
    )
}

combo_keys <- unique(scores_df[, c("lrt_subset", "method", "vector_label")])

for (i in seq_len(nrow(combo_keys))) {
  lrt_sub <- combo_keys$lrt_subset[i]
  method <- combo_keys$method[i]
  vec_label <- combo_keys$vector_label[i]

  sub_df <- scores_df[
    scores_df$lrt_subset == lrt_sub & scores_df$vector_label == vec_label,
  ]

  if (method == "best_fit") {
    title_str <- paste0("LRT (CPM+FC filtered): ", lrt_sub, " | Best-fit")
    fname <- paste0(lrt_sub, "_bestfit.svg")
  } else {
    title_str <- paste0(
      "LRT (CPM+FC filtered): ", lrt_sub, " | Vector: ",
      gsub("to", "->", vec_label)
    )
    fname <- paste0(lrt_sub, "_centroid_", vec_label, ".svg")
  }

  p <- make_panel(sub_df, title_str)

  svglite::svglite(file.path(fig_dir, fname),
    width = 4, height = 7,
    fix_text_size = FALSE
  )
  print(p)
  grDevices::dev.off()
}
