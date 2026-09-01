# ============================================================
# 01_fig2_sialo_deg_landscape.R
# Figure 2 -- sialylation-associated DEG landscape
#
# Pipeline:
#   (1) differential expression analysis (limma) across GEO cohorts
#   (2) intersection with the sialylation gene compendium
#   (3) univariate Cox screen of Sialo-DEGs in TCGA-HNSC
#   (4) GO/KEGG enrichment
#   (5) mutation landscape (oncoplot) and genomic distribution (circos)
#   (6) Figure 2 assembly (panels A-G)
#
# Requirements: see ENVIRONMENT.md
# ============================================================

rm(list = ls()); gc()

library(limma)
library(ggplot2)
library(pheatmap)
library(survival)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(maftools)
library(RCircos)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)

# palette used throughout the study
COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"

# ------------------------------------------------------------------
# (1) DEG analysis across GEO cohorts (batch-corrected)
# ------------------------------------------------------------------
# Raw microarray data were downloaded from GEO (accessions in Methods),
# normalised per platform, and probe-annotated; the processed matrices
# are stored under results/Rdata/. Only the modelling steps are shown.

cohort_files <- c(
  GSE29330 = "results/Rdata/GSE29330_preprocess.Rdata",
  GSE30784 = "results/Rdata/GSE30784_preprocess.Rdata",
  GSE3292  = "results/Rdata/GSE3292_preprocess.Rdata",
  GSE42743 = "results/Rdata/GSE42743_preprocess.Rdata",
  GSE6791  = "results/Rdata/GSE6791_preprocess.Rdata",
  GSE7224  = "results/Rdata/GSE7224_preprocess.Rdata",
  GSE9844  = "results/Rdata/GSE9844_preprocess.Rdata")

exps   <- vector("list", length(cohort_files))
groups <- vector("list", length(cohort_files))
for (i in seq_along(cohort_files)) {
  load(cohort_files[[i]])
  exps[[i]]   <- exprs_obj                     # normalised expression matrix
  groups[[i]] <- pdata$group                   # "Tumor" / "Normal" annotation
}

# restrict to genes measured in all cohorts, then merge
common_genes <- Reduce(intersect, lapply(exps, rownames))
expr_all <- do.call(cbind, lapply(exps, function(x) x[common_genes, , drop = FALSE]))
group  <- factor(unlist(groups), levels = c("Normal", "Tumor"))
batch  <- factor(rep(names(cohort_files), lengths(groups)))

# batch correction + empirical-Bayes moderated t-test
design <- model.matrix(~ group)
expr_corrected <- limma::removeBatchEffect(expr_all, batch = batch, design = design)
fit <- eBayes(lmFit(expr_corrected, design))
deg <- topTable(fit, coef = "groupTumor", number = Inf)
deg_sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ]

volcano_df <- deg
volcano_df$threshold <- ifelse(volcano_df$adj.P.Val < 0.05 & abs(volcano_df$logFC) > 1,
                               "Significant", "Not Significant")
ggplot(volcano_df, aes(logFC, -log10(adj.P.Val), color = threshold)) +
  geom_point(size = 1) + scale_color_manual(values = c("grey", "red")) +
  theme_bw() + ggtitle("Differential expression volcano")
ggsave("figures/Fig2_panelA_volcano.pdf", width = 6, height = 5)

# ------------------------------------------------------------------
# (2) intersection with the sialylation gene compendium
# ------------------------------------------------------------------

sialo_genes <- rownames(read.csv("data/sialylation_gene_compendium.csv", row.names = 1))
sialo_deg  <- intersect(rownames(deg_sig), sialo_genes)      # 156 Sialo-DEGs

venn <- ggvenn(list(`GEO DEGs` = rownames(deg_sig), `Sialylation genes` = sialo_genes),
               show_percentage = TRUE, fill_color = c(COL_NORMAL, COL_TUMOR))
ggsave("figures/Fig2_panelB_venn.pdf", venn, width = 5, height = 5)

# z-scored heatmap of the Sialo-DEGs (tumour vs normal)
expr_sub <- t(scale(t(expr_corrected[sialo_deg, ])))
expr_sub[expr_sub > 2] <- 2; expr_sub[expr_sub < -2] <- -2
ann_col <- data.frame(Group = group, row.names = colnames(expr_sub))
pheatmap(expr_sub, annotation_col = ann_col, show_colnames = FALSE,
         annotation_colors = list(Group = c(Normal = COL_NORMAL, Tumor = COL_TUMOR)),
         filename = "figures/Fig2_panelC_heatmap.pdf", height = 8)

# ------------------------------------------------------------------
# (3) univariate Cox screen of Sialo-DEGs in TCGA-HNSC
# ------------------------------------------------------------------

load("data/TCGA_HNSC_expression.Rdata")     # genes x samples TPM matrix
load("data/TCGA_HNSC_survival.Rdata")       # OS.time / OS
expr_tcga <- expr_tpm[sialo_deg, ]
expr_tcga <- expr_tcga[, rownames(clinical)]   # tumour samples with survival

cox_res <- do.call(rbind, lapply(rownames(expr_tcga), function(g) {
  fit <- coxph(Surv(OS.time, OS) ~ as.numeric(expr_tcga[g, ]), data = clinical)
  s <- summary(fit)$coefficients
  data.frame(gene = g, HR = exp(s[1]), pvalue = s[5])
}))
sig_genes <- cox_res$gene[cox_res$pvalue < 0.05]      # 31 prognostic Sialo-DEGs

# forest plot of the significant screen
cox_res$type <- ifelse(cox_res$HR > 1, "Risk", "Protect")
ggplot(cox_res[cox_res$pvalue < 0.05, ], aes(HR, reorder(gene, HR), color = type)) +
  geom_point(size = 2) + geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() + theme_bw() + labs(x = "Hazard ratio (log scale)", y = NULL)
ggsave("figures/Fig2_panelD_forest.pdf", width = 6, height = 8)

# ------------------------------------------------------------------
# (4) GO/KEGG enrichment of up-/down-regulated Sialo-DEGs
# ------------------------------------------------------------------

up_genes   <- sialo_deg[deg_sig[sialo_deg, "logFC"] > 0]
down_genes <- sialo_deg[deg_sig[sialo_deg, "logFC"] < 0]
up_entrez   <- bitr(up_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
down_entrez <- bitr(down_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID

ego_up   <- enrichGO(up_entrez, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
ego_down <- enrichGO(down_entrez, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
kegg_up   <- enrichKEGG(up_entrez, organism = "hsa")
kegg_down <- enrichKEGG(down_entrez, organism = "hsa")

enrich_df <- bind_rows(
  as.data.frame(ego_up)   %>% mutate(Regulation = "Up",   Category = "GO"),
  as.data.frame(ego_down) %>% mutate(Regulation = "Down", Category = "GO"),
  as.data.frame(kegg_up)  %>% mutate(Regulation = "Up",   Category = "KEGG"),
  as.data.frame(kegg_down) %>% mutate(Regulation = "Down", Category = "KEGG"))
top_enrich <- enrich_df %>% group_by(Category, Regulation) %>%
  arrange(pvalue) %>% slice_head(n = 8)

ggplot(top_enrich, aes(-log10(pvalue), reorder(Description, -log10(pvalue)), fill = Regulation)) +
  geom_col(position = "dodge") + facet_grid(Category ~ ., scales = "free_y") +
  theme_bw() + labs(x = "-log10(p-value)", y = NULL)
ggsave("figures/Fig2_panelE_enrichment.pdf", width = 8, height = 8)

# ------------------------------------------------------------------
# (5) mutation landscape and genomic distribution
# ------------------------------------------------------------------

# MAF obtained from the TCGA MC3 project; restricted to the 31 genes
maf <- read.delim(gzfile("data/TCGA_HNSC_mc3.maf.gz"), comment.char = "#")
maf <- maf[maf$Hugo_Symbol %in% sig_genes, ]
maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 12)

pdf("figures/Fig2_panelG_oncoplot.pdf", width = 8, height = 6)
maftools::oncoplot(read.maf(maf), genes = sig_genes,
                   showTumorSampleBarcodes = FALSE, removeNonMutated = TRUE)
dev.off()

# circos: genomic coordinates retrieved via biomaRt (not shown; cached in
# results/plot_data/), plotted as gene tiles around the ideogram
load("results/plot_data/gene_circos.Rdata")
pdf("figures/Fig2_panelF_circos.pdf", width = 6, height = 6)
RCircos.Set.Core.Components(cyto_info, tracks.inside = 3, tracks.outside = 0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Tile.Plot(chr_ring, track.num = 1, side = "in")
RCircos.Gene.Connector.Plot(label_data, track.num = 2, side = "in")
RCircos.Gene.Name.Plot(label_data, name.col = 4, track.num = 3, side = "in")
text(0, 0, "Homo sapiens", cex = 1.4)
dev.off()

# ------------------------------------------------------------------
# (6) Figure 2 assembly
# ------------------------------------------------------------------
# Panels A-G are combined from the PDFs above with the layout defined in
# the manuscript; representative layout code:
#   library(patchwork)
#   (panelA | panelB) / (panelC | panelD) / (panelE | panelF | panelG)
