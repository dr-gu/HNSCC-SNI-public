# ============================================================
# 02_fig3_molecular_subtypes.R
# Figure 3 -- sialylation-defined molecular subtypes (C1/C2)
#
# Pipeline:
#   (1) consensus clustering of the prognostic Sialo-DEGs (k = 2)
#   (2) Kaplan-Meier survival of the two clusters
#   (3) GSEA between clusters (KEGG + GO)
#   (4) ESTIMATE scores (tumour purity / immune / stromal)
#   (5) immune-checkpoint gene expression
#   (6) TIP cancer-immunity cycle steps (ssGSEA)
#   (7) Figure 3 assembly (panels A-F)
#
# Requirements: see ENVIRONMENT.md
# ============================================================

rm(list = ls()); gc()

library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(estimate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"

# ------------------------------------------------------------------
# (1) consensus clustering of the 31 prognostic Sialo-DEGs
# ------------------------------------------------------------------
# Input matrix: TCGA-HNSC expression of the genes retained by the
# univariate Cox screen (script 01), rows = genes, columns = samples.

expr <- read.delim("data/sialo_prognostic_expression.txt", row.names = 1, check.names = FALSE)
expr_scaled <- t(scale(t(as.matrix(expr))))

cc <- ConsensusClusterPlus(
  expr_scaled, maxK = 10, reps = 1000, pItem = 0.8, pFeature = 1,
  clusterAlg = "hc", distance = "pearson", seed = 42,
  plot = "pdf", title = "figures/consensus_cluster")
save(cc, file = "results/consensus_cluster.Rdata")

cluster <- cc[[2]]$consensusClass          # k = 2: C1 / C2
write.table(data.frame(sample = names(cluster), cluster = cluster),
            "tables/cluster_assignment.txt", sep = "\t", row.names = FALSE)

# consensus-matrix heatmap (Fig. 3a)
m <- cc[[2]]$consensusMatrix
m <- m[cc[[2]]$consensusTree$order, cc[[2]]$consensusTree$order]
ann <- data.frame(Cluster = paste0("C", cluster[colnames(m)]), row.names = colnames(m))
pdf("figures/Fig3_panelA_consensus.pdf", width = 5, height = 4)
Heatmap(m, col = circlize::colorRamp2(c(0, 1), c("white", "#3B528B")),
        cluster_rows = FALSE, cluster_columns = FALSE,
        top_annotation = HeatmapAnnotation(Cluster = ann$Cluster,
          col = list(Cluster = c(C1 = "#FD8D3C", C2 = "#74ADD1"))),
        show_row_names = FALSE, show_column_names = FALSE)
dev.off()

# ------------------------------------------------------------------
# (2) Kaplan-Meier of the two clusters (Fig. 3b)
# ------------------------------------------------------------------

clinical <- read.delim("data/tcga_hnsc_clinical.txt", row.names = 1)   # OS.time / OS
clinical$cluster <- factor(cluster[rownames(clinical)], levels = 1:2, labels = c("C1", "C2"))

fit <- survfit(Surv(OS.time, OS) ~ cluster, data = clinical)
ggsurvplot(fit, data = clinical, pval = TRUE, risk.table = TRUE,
           palette = c(COL_TUMOR, COL_NORMAL), legend.labs = c("C1", "C2"),
           filename = "figures/Fig3_panelB_km.pdf")

# ------------------------------------------------------------------
# (3) GSEA between clusters (Fig. 3c)
# ------------------------------------------------------------------

load("data/TCGA_HNSC_expression.Rdata")   # genes x samples TPM matrix
expr_tcga <- expr_tpm[, rownames(clinical)]
design <- model.matrix(~ cluster)
fit <- eBayes(lmFit(expr_tcga, design))
gene_list <- sort(setNames(fit$t[, 2], rownames(fit$t)), decreasing = TRUE)
entrez <- bitr(names(gene_list), "SYMBOL", "ENTREZID", org.Hs.eg.db)
gene_list <- sort(setNames(gene_list[entrez$SYMBOL], entrez$ENTREZID), decreasing = TRUE)

gsea_kegg <- gseKEGG(geneList = gene_list, organism = "hsa", pvalueCutoff = 0.05)
gsea_go   <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP",
                   pvalueCutoff = 0.05)

top_sets <- as.data.frame(gsea_kegg) %>% filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>% slice_head(n = 20)
ggplot(top_sets, aes(NES, reorder(Description, NES), size = setSize, color = p.adjust)) +
  geom_point() + scale_color_gradient(low = COL_TUMOR, high = COL_NORMAL) +
  theme_bw() + labs(x = "Normalised enrichment score", y = NULL)
ggsave("figures/Fig3_panelC_gsea.pdf", width = 7, height = 5)

# ------------------------------------------------------------------
# (4) ESTIMATE scores (Fig. 3d)
# ------------------------------------------------------------------

write.table(expr_tcga, "estimate_input.gct", sep = "\t", quote = FALSE)
filterCommonGenes("estimate_input.gct", "estimate_genes.gct", id = "GeneSymbol")
estimateScore("estimate_genes.gct", "estimate_score.gct")
est <- read.table("estimate_score.gct", skip = 2, header = TRUE, row.names = 1)
est <- as.data.frame(t(est[, -1]))
est$cluster <- cluster[rownames(est)]

est_long <- est %>% pivot_longer(c(StromalScore, ImmuneScore, ESTIMATEScore),
                                 names_to = "Score", values_to = "value")
ggboxplot(est_long, x = "cluster", y = "value", fill = "cluster",
          palette = c(COL_TUMOR, COL_NORMAL), add = "jitter", facet.by = "Score",
          scales = "free") + stat_compare_means(method = "wilcox.test") +
  labs(x = NULL, y = "Score")
ggsave("figures/Fig3_panelD_estimate.pdf", width = 9, height = 4)

# ------------------------------------------------------------------
# (5) immune-checkpoint gene expression (Fig. 3e)
# ------------------------------------------------------------------

checkpoints <- c("CD274", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "PDCD1LG2", "SIGLEC15", "TIGIT")
ckp <- as.data.frame(t(log2(expr_tcga[checkpoints, ] + 1)))
ckp$cluster <- cluster[rownames(ckp)]
ckp_long <- ckp %>% pivot_longer(all_of(checkpoints), names_to = "gene", values_to = "expr")
ggboxplot(ckp_long, x = "cluster", y = "expr", fill = "cluster",
          palette = c(COL_TUMOR, COL_NORMAL), facet.by = "gene", nrow = 1,
          scales = "free_y") +
  stat_compare_means(comparisons = list(c("C1", "C2")), method = "wilcox.test",
                     label = "p.signif", bracket.size = 0) +
  labs(x = NULL, y = "log2(TPM+1)")
ggsave("figures/Fig3_panelE_checkpoints.pdf", width = 12, height = 4)

# ------------------------------------------------------------------
# (6) TIP cancer-immunity cycle steps (Fig. 3f)
# ------------------------------------------------------------------
# TIP ssGSEA scoring is performed by the third-party TIP tool
# (see tools/); its per-sample normalised scores are loaded here.

load("results/tip_ssgsea_scores.Rdata")      # steps x samples matrix
tip <- as.data.frame(t(tip_scores))
tip$cluster <- cluster[rownames(tip)]
tip_long <- tip %>% pivot_longer(starts_with("Step"), names_to = "Step", values_to = "Score")
tip_long$Score <- pmin(pmax(scale(tip_long$Score), -3), 3)

ggboxplot(tip_long, x = "cluster", y = "Score", fill = "cluster",
          palette = c(COL_TUMOR, COL_NORMAL), facet.by = "Step", nrow = 1) +
  stat_compare_means(comparisons = list(c("C1", "C2")), method = "wilcox.test",
                     label = "p.signif", bracket.size = 0) +
  labs(x = NULL, y = "Scaled ssGSEA score")
ggsave("figures/Fig3_panelF_tip.pdf", width = 14, height = 4)

# ------------------------------------------------------------------
# (7) Figure 3 assembly
# ------------------------------------------------------------------
# Panels A-F are combined from the PDFs above with the layout defined in
# the manuscript; representative layout code:
#   library(patchwork)
#   (panelA | panelB | panelC) / (panelD / panelE / panelF)
