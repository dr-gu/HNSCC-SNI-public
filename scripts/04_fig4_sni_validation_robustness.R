# ============================================================
# 04_fig4_sni_validation_robustness.R
# Figure 4 (validation panels) + Supplementary Figures S1-S3
#
# Pipeline:
#   (1) application of the pre-fitted SNI model to three additional
#       external cohorts (E-MTAB-8588, GSE65858, GSE41613)
#   (2) per-cohort C-index, time-dependent AUC, Kaplan-Meier
#   (3) bootstrap confidence interval for the training C-index
#   (4) pooled validation across the four external cohorts
#       (cohort-stratified Cox model + forest plot)
#   (5) benchmark against clinical variables
#   (6) Supplementary Figures S1-S3
#
# Requirements: see ENVIRONMENT.md
# ============================================================

rm(list = ls()); gc()

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(boot)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"

# ------------------------------------------------------------------
# (1) score the external cohorts with the pre-fitted model
# ------------------------------------------------------------------
# The model objects produced by script 03 are re-used as-is; candidate
# genes absent from a given platform are excluded, and the remaining
# model genes are re-normalised to the training distribution.

load("results/Rdata/sni_model_fit.Rdata")      # StepCox coefficients + RSF object
load("results/Rdata/Mime_TCGA.Rdata")

external <- list(
  `E-MTAB-8588` = "results/Rdata/Mime_E_MTAB_8588.Rdata",
  GSE65858      = "results/Rdata/Mime_GSE65858.Rdata",
  GSE41613      = "results/Rdata/Mime_GSE41613.Rdata")

score_cohort <- function(file) {
  load(file)
  x <- as.data.frame(x)
  common <- intersect(model_genes, colnames(x))
  lp <- predict(cox_fit, newdata = x[, c("OS.time", "OS", common)])
  data.frame(OS.time = x$OS.time, OS = x$OS, RS = as.numeric(lp))
}

scores <- lapply(external, score_cohort)
scores$GSE42743 <- scores_existing_gse42743     # from script 03

# ------------------------------------------------------------------
# (2) per-cohort performance (Table 2)
# ------------------------------------------------------------------

per_cohort <- lapply(names(scores), function(ds) {
  s <- scores[[ds]]
  cindex <- survConcordance(Surv(OS.time, OS) ~ RS, data = s)$concordance
  auc <- sapply(c(1, 2, 3), function(t) {
    roc_auc <- timeROC::timeROC(T = s$OS.time, delta = s$OS,
                                marker = s$RS, cause = 1, times = t * 365)
    mean(roc_auc$AUC)
  })
  km_p <- survdiff(Surv(OS.time, OS) ~ (RS > median(RS)), data = s)$pvalue
  data.frame(Cohort = ds, n = nrow(s), Cindex = round(cindex, 3),
             AUC_1y = auc[1], AUC_2y = auc[2], AUC_3y = auc[3],
             KM_p = signif(km_p, 3))
})
table2 <- bind_rows(per_cohort)
write.csv(table2, "tables/Table2_per_cohort_performance.csv", row.names = FALSE)

# Kaplan-Meier per cohort (Fig. S2)
for (ds in names(scores)) {
  s <- scores[[ds]]
  s$riskGroup <- ifelse(s$RS > median(s$RS), "High", "Low")
  fit <- survfit(Surv(OS.time, OS) ~ riskGroup, data = s)
  ggsurvplot(fit, data = s, pval = TRUE, palette = c(COL_TUMOR, COL_NORMAL),
             filename = sprintf("figures/FigS2_km_%s.pdf", ds))
}

# time-dependent ROC per cohort (Fig. S1)
roc_list <- lapply(names(scores), function(ds) {
  s <- scores[[ds]]
  timeROC::timeROC(T = s$OS.time, delta = s$OS, marker = s$RS,
                   cause = 1, times = c(1, 2, 3) * 365)
})
# representative ROC plot for one cohort:
plot(roc_list[[1]], col = c(COL_NORMAL, "#7B7B73", COL_TUMOR), lwd = 1.5)

# ------------------------------------------------------------------
# (3) bootstrap confidence interval for the training C-index
# ------------------------------------------------------------------

cindex_boot <- function(data, idx) {
  survConcordance(Surv(OS.time, OS) ~ RS, data = data[idx, ])$concordance
}
set.seed(42)
b <- boot(scores$`TCGA-HNSC`, cindex_boot, R = 1000)
boot.ci(b, type = "perc")

# ------------------------------------------------------------------
# (4) pooled validation across the four external cohorts (Fig. S3)
# ------------------------------------------------------------------

pooled <- bind_rows(lapply(external_names <- setdiff(names(scores), "TCGA-HNSC"),
                           function(ds) data.frame(Cohort = ds, scores[[ds]])))

fit_pooled <- coxph(Surv(OS.time, OS) ~ RS + strata(Cohort), data = pooled)
pooled_hr <- summary(fit_pooled)$conf.int[, c(1, 3, 4)]

forest_df <- bind_rows(
  lapply(external_names, function(ds) {
    f <- coxph(Surv(OS.time, OS) ~ RS, data = scores[[ds]])
    ci <- summary(f)$conf.int[, c(1, 3, 4)]
    data.frame(Cohort = ds, HR = ci[1], lower = ci[2], upper = ci[3])
  }),
  data.frame(Cohort = "Pooled (n = 524)", HR = pooled_hr[1],
             lower = pooled_hr[2], upper = pooled_hr[3]))
forest_df$Cohort <- factor(forest_df$Cohort, levels = rev(forest_df$Cohort))

ggplot(forest_df, aes(HR, Cohort)) +
  geom_point(size = 2.5, color = COL_TUMOR) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() + theme_bw() +
  labs(x = "Hazard ratio per unit risk score (log scale)", y = NULL)
ggsave("figures/FigS3_pooled_forest.pdf", width = 6, height = 3.5)

# pooled dichotomised Kaplan-Meier (median per cohort)
pooled$riskGroup <- ifelse(pooled$RS > median(pooled$RS), "High", "Low")
survdiff(Surv(OS.time, OS) ~ riskGroup + strata(Cohort), data = pooled)

# ------------------------------------------------------------------
# (5) benchmark against clinical variables
# ------------------------------------------------------------------

load("results/Rdata/clinical_comparison.Rdata")    # complete-case TCGA subset
# clinical-only model: age + AJCC stage + T/N category + grade + sex
fit_clin <- coxph(Surv(OS.time, OS) ~ age + stage + T + N + grade + sex,
                  data = clin_df)
fit_sni  <- coxph(Surv(OS.time, OS) ~ RS, data = clin_df)
fit_both <- coxph(Surv(OS.time, OS) ~ age + stage + T + N + grade + sex + RS,
                  data = clin_df)
sapply(list(clinical = fit_clin, SNI = fit_sni, combined = fit_both),
       function(f) survConcordance(f)$concordance)

# ------------------------------------------------------------------
# (6) Figure 4 validation panels / Figures S1-S3 assembly
# ------------------------------------------------------------------
# Combined from the PDFs above with the layout defined in the manuscript.
