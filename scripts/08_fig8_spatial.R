# ============================================================
# 08_fig8_spatial.R
# Figure 8 -- spatial transcriptomics (10x Visium)
#
# Pipeline:
#   (1) spatial data loading and processing (Seurat)
#   (2) spatial distribution of the key genes and niche scores
#   (3) Figure 8 assembly
#
# Requirements: see ENVIRONMENT.md (spatial analysis uses the Seurat
# version noted there)
# ============================================================

rm(list = ls()); gc()

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"
KEY_GENES  <- c("HSPH1", "ST6GALNAC1")

# ------------------------------------------------------------------
# (1) load and process 10x Visium data
# ------------------------------------------------------------------
# Space Ranger output (filtered_feature_bc_matrix + spatial/) for each
# HNSCC section; sections are processed independently and then merged.

section_dirs <- list.files("data/visium/", pattern = "^S\\d+$", full.names = TRUE)

load_section <- function(d) {
  s <- Load10X_Spatial(data.dir = file.path(d, "outs"))
  s$section <- basename(d)
  s <- SCTransform(s, assay = "Spatial", verbose = FALSE)
  s
}

sections <- lapply(section_dirs, load_section)
st <- merge(sections[[1]], sections[-1])

# ------------------------------------------------------------------
# (2) spatial feature plots (Fig. 8a/b)
# ------------------------------------------------------------------

for (g in KEY_GENES) {
  SpatialFeaturePlot(st, features = g, alpha = c(0.1, 1), ncol = 2)
  ggsave(sprintf("figures/Fig8_feature_%s.pdf", g), width = 12, height = 8)
}

# niche score: mean z-scored expression of the model genes per spot
model_genes <- readLines("data/model_genes.txt")
st <- AddModuleScore(st, features = list(model_genes), name = "SNI_niche")
SpatialFeaturePlot(st, features = "SNI_niche1", alpha = c(0.1, 1))
ggsave("figures/Fig8_panelC_niche_score.pdf", width = 10, height = 7)

# ------------------------------------------------------------------
# (3) spot-level colocalisation of the two key genes
# ------------------------------------------------------------------
# Representative correlation of the two key genes across spots of one
# section, illustrating their co-enrichment within the tumour niche.

sec <- sections[[1]]
expr <- GetAssayData(sec, slot = "data")[KEY_GENES, ]
plot_df <- data.frame(HSPH1 = expr["HSPH1", ], ST6GALNAC1 = expr["ST6GALNAC1", ])
ggplot(plot_df, aes(HSPH1, ST6GALNAC1)) +
  geom_point(size = 0.8, alpha = 0.6, color = COL_TUMOR) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30") +
  ggpubr::stat_cor(method = "spearman") +
  theme_bw() + labs(x = "HSPH1", y = "ST6GALNAC1")
ggsave("figures/Fig8_panelD_colocalisation.pdf", width = 5, height = 4)

# ------------------------------------------------------------------
# (4) Figure 8 assembly
# ------------------------------------------------------------------
# Panels A-D are combined from the PDFs above with the layout defined
# in the manuscript.
