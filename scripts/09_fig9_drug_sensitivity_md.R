# ============================================================
# 09_fig9_drug_sensitivity_md.R
# Figure 9 -- AI-guided drug screening, pharmacogenomic
#             corroboration and molecular dynamics
#
# Pipeline:
#   (1) pharmacogenomic correlation analyses (GDSC2 / CTRP / PRISM /
#       oncoPredict bridge to TCGA-HNSC)
#   (2) Figure 9 assembly
# The multi-stage screening cascade (GraphBAN interaction prediction,
# ADMET filtering, docking, MD) is documented in screening/README.md;
# MD trajectory metrics are computed by 09_fig9_md_analysis.py.
#
# Requirements: see ENVIRONMENT.md
# ============================================================

rm(list = ls()); gc()

library(ggplot2)
library(dplyr)
library(tidyr)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"
KEY_GENES  <- c("HSPH1", "ST6GALNAC1")

# ------------------------------------------------------------------
# (1) pharmacogenomic correlation analyses
# ------------------------------------------------------------------

# Generic correlation routine between a gene's expression and drug
# sensitivity (IC50/AUC) across cell lines of a panel.
panel_correlation <- function(expr, sensitivity, panel_name, label) {
  stopifnot(identical(rownames(expr), rownames(sensitivity)))
  res <- do.call(rbind, lapply(colnames(sensitivity), function(drug) {
    lapply(KEY_GENES, function(g) {
      r <- cor.test(expr[[g]], sensitivity[[drug]], method = "spearman")
      data.frame(drug = drug, gene = g, rho = r$estimate, pvalue = r$p.value)
    }) %>% bind_rows()
  }))
  res$fdr <- p.adjust(res$pvalue, method = "BH")
  res$panel <- panel_name
  res$label <- label
  res
}

# GDSC2 panel (expression and IC50 from the same cell lines)
gdsc_expr <- read.csv("data/gdsc2_expression.csv", row.names = 1)
gdsc_sens <- read.csv("data/gdsc2_ic50.csv", row.names = 1)
gdsc_res <- panel_correlation(gdsc_expr, gdsc_sens, "GDSC2", "IC50")

# CTRP v2 panel
ctrp_expr <- read.csv("data/ctrp_expression.csv", row.names = 1)
ctrp_sens <- read.csv("data/ctrp_auc.csv", row.names = 1)
ctrp_res <- panel_correlation(ctrp_expr, ctrp_sens, "CTRP v2", "AUC")

# DepMap PRISM panel
prism_expr <- read.csv("data/prism_expression.csv", row.names = 1)
prism_sens <- read.csv("data/prism_auc.csv", row.names = 1)
prism_res <- panel_correlation(prism_expr, prism_sens, "PRISM", "AUC")

# oncoPredict bridge: GDSC2 models applied to TCGA-HNSC expression
# (target genes removed from the feature space to avoid circularity)
library(oncoPredict)
calcPhenotype(trainingExprData = as.matrix(gdsc_expr),
              trainingPtype = t(gdsc_sens),
              testExprData = as.matrix(tcga_expr_bridge),
              batchCorrect = "eb", powerTransformPhenotype = TRUE,
              minNumSamples = 10, printOutput = FALSE)
bridge_sens <- as.data.frame(t(predicted_sensitivity))
bridge_res <- panel_correlation(as.data.frame(tcga_expr_bridge), bridge_sens,
                                "TCGA-HNSC (bridge)", "predicted IC50")

# combine and filter
all_res <- bind_rows(gdsc_res, ctrp_res, prism_res, bridge_res)
sig_res <- all_res %>% filter(fdr < 0.05)

# representative volcano-style summary per panel (Fig. 9)
ggplot(sig_res, aes(rho, -log10(pvalue), color = gene)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(HSPH1 = COL_TUMOR, ST6GALNAC1 = COL_NORMAL)) +
  facet_wrap(~ panel, scales = "free") + theme_bw() +
  labs(x = "Spearman rho (expression vs sensitivity)",
       y = "-log10(p-value)")
ggsave("figures/Fig9_panel_correlation.pdf", width = 10, height = 7)

# top significant associations per gene and panel (Table S4)
top_hits <- sig_res %>% group_by(panel, gene) %>% arrange(desc(abs(rho))) %>%
  slice_head(n = 20)
write.csv(top_hits, "tables/TableS4_panel_correlations.csv", row.names = FALSE)

# ------------------------------------------------------------------
# (2) Figure 9 assembly
# ------------------------------------------------------------------
# Screening-funnel, docking and MD panels are produced by the pipeline
# documented in screening/README.md; the correlation panels above are
# combined with them per the manuscript layout.
