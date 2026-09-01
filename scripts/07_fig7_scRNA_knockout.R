# ============================================================
# 07_fig7_scRNA_knockout.R
# Figure 7 -- single-cell atlas and virtual knockout of the key genes
#
# Pipeline:
#   (1) single-cell processing (Seurat standard workflow)
#   (2) cell-type annotation and key-gene distribution
#   (3) virtual knockout of HSPH1 / ST6GALNAC1 (scTenifoldKnk)
#   (4) Figure 7 assembly
#
# Requirements: see ENVIRONMENT.md (scRNA analysis uses the Seurat
# version noted there)
# ============================================================

rm(list = ls()); gc()

library(Seurat)
library(ggplot2)
library(dplyr)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"
KEY_GENES  <- c("HSPH1", "ST6GALNAC1")

# ------------------------------------------------------------------
# (1) single-cell processing
# ------------------------------------------------------------------
# Pre-annotated Seurat objects of the HNSCC scRNA-seq dataset were
# generated externally; the workflow below reproduces the standard
# quality control, normalisation and clustering steps.

sc <- readRDS("data/scRNA_seurat_object.rds")    # already QC-filtered

sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
sc <- ScaleData(sc, features = VariableFeatures(sc))
sc <- RunPCA(sc, npcs = 30)
sc <- FindNeighbors(sc, dims = 1:30)
sc <- FindClusters(sc, resolution = 0.5)
sc <- RunUMAP(sc, dims = 1:30)

DimPlot(sc, group.by = "celltype", label = TRUE, repel = TRUE)
ggsave("figures/Fig7_panelA_umap.pdf", width = 7, height = 6)

# ------------------------------------------------------------------
# (2) key-gene distribution across cell types (Fig. 7b)
# ------------------------------------------------------------------

FeaturePlot(sc, features = KEY_GENES, split.by = "celltype", ncol = 2)
ggsave("figures/Fig7_panelB_feature.pdf", width = 10, height = 8)

VlnPlot(sc, features = KEY_GENES, group.by = "celltype",
        pt.size = 0, cols = scales::hue_pal()(length(unique(sc$celltype))))
ggsave("figures/Fig7_panelC_vln.pdf", width = 12, height = 5)

# ------------------------------------------------------------------
# (3) virtual knockout of the key genes (scTenifoldKnk)
# ------------------------------------------------------------------
# The scTenifoldKnk framework compares the observed gene-regulatory
# network against a network with the target gene virtually knocked
# out, and ranks the most affected genes.

library(scTenifoldKnk)

counts <- GetAssayData(sc, slot = "counts")

knk_hsph1 <- scTenifoldKnk(counts, gKO = "HSPH1", qc_minLibSize = 100,
                           qc_minCells = 10)
knk_st6   <- scTenifoldKnk(counts, gKO = "ST6GALNAC1", qc_minLibSize = 100,
                           qc_minCells = 10)

# distance-ranked gene list for each knockout
rank_genes <- function(knk, n = 50) {
  df <- knk$diffRegulation
  df <- df[order(df$distance, decreasing = TRUE), ]
  head(df$gene, n)
}
top_hsph1 <- rank_genes(knk_hsph1)
top_st6   <- rank_genes(knk_st6)

# representative visualisation of the KO-affected network (Fig. 7d)
plotKO <- function(knk, genes) {
  g <- knk$tensorNetworks[[1]]
  sub <- igraph::induced_subgraph(g, igraph::V(g)[igraph::V(g)$name %in% genes])
  plot(sub, vertex.size = 4, vertex.label.cex = 0.6,
       vertex.color = ifelse(igraph::V(sub)$name %in% KEY_GENES, COL_TUMOR, "grey80"))
}
pdf("figures/Fig7_panelD_knockout_network.pdf", width = 8, height = 8)
plotKO(knk_hsph1, top_hsph1)
dev.off()

save(knk_hsph1, knk_st6, file = "results/virtual_knockout.Rdata")

# ------------------------------------------------------------------
# (4) Figure 7 assembly
# ------------------------------------------------------------------
# Panels A-D are combined from the PDFs above with the layout defined
# in the manuscript.
