# ============================================================
# 03_fig4_sni_construction.R
# Figure 4 (construction panels) -- Sialylation Niche Index (SNI)
#
# Pipeline:
#   (1) training/validation data preparation
#   (2) benchmark of 101 machine-learning algorithm combinations (Mime1)
#   (3) model selection by the mean C-index across the two cohorts
#   (4) risk-score derivation, PCA separation, Kaplan-Meier,
#       time-dependent ROC, risk-factor scatter plots
#   (5) Figure 4 construction panels
#
# Requirements: see ENVIRONMENT.md
# ============================================================

rm(list = ls()); gc()

library(Mime1)
library(survival)
library(ggplot2)
library(patchwork)
library(scatterplot3d)
library(dplyr)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"

# ------------------------------------------------------------------
# (1) data preparation
# ------------------------------------------------------------------
# Training cohort: TCGA-HNSC (RNA-seq, n = 501 with overall survival).
# Validation cohort: GSE42743 (microarray, n = 74 with overall survival).
# Raw data were downloaded, normalised, probe-collapsed and matched to
# the candidate genes as described in Methods; the ready matrices are
# stored under results/Rdata/.

load("results/Rdata/Mime_TCGA.Rdata")        # futime / fustat / genes
load("results/Rdata/Mime_GSE42743.Rdata")

datasets <- list(`TCGA-HNSC` = Mime_TCGA, GSE42743 = GSE42743_survival)
datasets <- lapply(datasets, function(x) {
  x <- as.data.frame(x)
  x <- cbind(ID = rownames(x), x)
  colnames(x)[1:3] <- c("ID", "OS.time", "OS")
  x[, -c(1:3)] <- scale(x[, -c(1:3)])
  x
})

candidate_genes <- read.delim("data/prognostic_sialo_genes.txt", header = FALSE)[, 1]

# ------------------------------------------------------------------
# (2) benchmark of 101 algorithm combinations
# ------------------------------------------------------------------
# The Mime1 framework benchmarks combinations of 10 machine-learning
# algorithms (survival-specific) under a 3-fold cross-validated tuning
# scheme, yielding 101 unique combinations; the development call is:

res <- ML.Dev.Prog.Sig(
  train_data = datasets$`TCGA-HNSC`,
  list_train_vali_Data = datasets,
  unicox.filter.for.candi = TRUE,
  unicox_p_cutoff = 0.01,
  candidate_genes = candidate_genes,
  mode = "all",            # run all 101 algorithm combinations
  nodesize = 10,
  seed = 42)

save(res, file = "results/Rdata/model_benchmark.Rdata")

# per-combination C-index overview (Fig. 4a)
rank_df <- res$Cindex.res %>%
  mutate(mean_cindex = rowMeans(across(starts_with("Cindex")))) %>%
  arrange(desc(mean_cindex))
ggplot(rank_df, aes(mean_cindex, reorder(Model, mean_cindex))) +
  geom_point() + theme_bw() +
  labs(x = "Mean C-index (training + validation)", y = NULL) +
  geom_point(data = rank_df[1, ], color = COL_TUMOR, size = 3)
ggsave("figures/Fig4_panelA_ranking.pdf", width = 5, height = 12)

# ------------------------------------------------------------------
# (3) model selection
# ------------------------------------------------------------------
# Selection criterion: highest mean C-index across the training and
# validation cohorts; the selected combination was used thereafter
# without modification.

best_model <- rank_df$Model[1]                      # StepCox[forward] + RSF
risk_all <- res$riskscore[[best_model]]

# ------------------------------------------------------------------
# (4) risk score, survival and discrimination (Fig. 4b-j)
# ------------------------------------------------------------------

for (ds in names(risk_all)) {
  risk_df <- risk_all[[ds]]
  risk_df$riskGroup <- ifelse(risk_df$RS > median(risk_df$RS), "High", "Low")

  # PCA separation of the two risk groups (Fig. 4b/c)
  expr_sub <- datasets[[ds]][, res$Sig.genes]
  pca <- prcomp(expr_sub, scale. = TRUE)
  pca_df <- data.frame(pca$x[, 1:3], group = risk_df$riskGroup)
  pdf(sprintf("figures/Fig4_pca_%s.pdf", ds), width = 6, height = 6)
  scatterplot3d(pca_df[, 1:3], color = ifelse(pca_df$group == "High", COL_TUMOR, COL_NORMAL),
                pch = 16, main = sprintf("3D PCA of %s", ds))
  dev.off()

  # Kaplan-Meier of the risk groups (Fig. 4d/e)
  fit <- survfit(Surv(OS.time, OS) ~ riskGroup, data = risk_df)
  p_km <- survminer::ggsurvplot(fit, data = risk_df, pval = TRUE,
                                palette = c(COL_TUMOR, COL_NORMAL))$plot
  ggsave(sprintf("figures/Fig4_km_%s.pdf", ds), p_km, width = 5, height = 4)

  # time-dependent ROC at 1/2/3 years (Fig. 4f/g)
  aucs <- lapply(c(1, 2, 3), function(t) {
    cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,
                   train_data = datasets$`TCGA-HNSC`,
                   inputmatrix.list = datasets,
                   mode = "all", AUC_time = t, auc_cal_method = "KM")
  })
  roc_vis(aucs[[3]], model_name = best_model, dataset = ds,
          order = ds, anno_position = c(0.65, 0.55), year = 3)

  # risk-score distribution and survival tracking (Fig. 4h-j)
  risk_df <- risk_df[order(risk_df$RS), ]
  risk_df$index <- seq_len(nrow(risk_df))
  p1 <- ggplot(risk_df, aes(index, RS, color = riskGroup)) +
    geom_point(size = 1.2) + scale_color_manual(values = c(COL_TUMOR, COL_NORMAL)) +
    theme_classic() + labs(y = "Risk score", x = NULL)
  p2 <- ggplot(risk_df, aes(index, OS.time, color = factor(OS))) +
    geom_point(size = 1.2) + scale_color_manual(values = c(COL_NORMAL, COL_TUMOR)) +
    theme_classic() + labs(y = "Survival time (days)", x = NULL)
  ggsave(sprintf("figures/Fig4_tracking_%s.pdf", ds), p1 / p2, width = 4, height = 7)
}

save(risk_all, res$Sig.genes, file = "results/Rdata/sni_risk_scores.Rdata")

# ------------------------------------------------------------------
# (5) Figure 4 assembly (construction panels)
# ------------------------------------------------------------------
# The construction panels are combined from the PDFs above; the external
# validation panels are produced by script 04.
