# Differentiation

PCA-based method for quantifying differentiation
timing differences between conditions in bulk RNA-seq
time course data.

**[View the full interactive analysis](https://zohebkhan1.github.io/differentiation/)**

## Method

1. Identify temporally variable genes in a control
   time course via DESeq2 likelihood ratio test
2. Build a PCA space from the top LRT genes using
   control samples only
3. Define a reference vector between early and late
   control centroids in PC space
4. Project all samples onto this vector to compute
   a normalized maturation score (0 = early, 1 = late)

![Method overview](assets/overview.svg)

## Dataset

Publicly available data from
[Martinez et al. (2024)](https://doi.org/10.3389/fncel.2024.1341141).
18 samples: D21 (euploid control) and T21 (trisomy 21)
iPSC-derived neural progenitors at Days 6, 10, and 17.

## Repository Structure

```
dat/              Input data (counts, metadata)
assets/           SVG diagrams
lrt/              Reference LRT scripts (DSHC1 project)
index.qmd         Quarto analysis notebook
_quarto.yml       Quarto project config
docs/             Rendered GitHub Pages site
```

## Dependencies

- R (>= 4.0)
- DESeq2, ComplexHeatmap, plotly, ggplot2, dplyr,
  tidyr, tibble, circlize
- Quarto (>= 1.4) for rendering

## Rendering

```bash
quarto render
```

Output is written to `docs/` and served via
GitHub Pages.
