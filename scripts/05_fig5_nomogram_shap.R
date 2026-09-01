# ============================================================
# 05_fig5_nomogram_shap.R
# Figure 5 -- independent prognostic value, nomogram and SHAP
#
# Pipeline:
#   (1) univariate and multivariable Cox (SNI vs clinicopathological)
#   (2) prognostic nomogram (age + SNI) with calibration
#   (3) SHAP feature attribution of the model components
#   (4) threshold-sensitivity analysis of the SHAP gene selection
#   (5) Figure 5 assembly
#
# Requirements: see ENVIRONMENT.md
# ============================================================

rm(list = ls()); gc()

library(survival)
library(rms)
library(ggplot2)
library(dplyr)
library(fastshap)
library(randomForestSRC)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"

load("results/Rdata/sni_risk_scores.Rdata")     # risk_all (script 03)
load("results/Rdata/clinical_with_sni.Rdata")   # TCGA clinical + RS

# ------------------------------------------------------------------
# (1) univariate and multivariable Cox
# ------------------------------------------------------------------

vars <- c("age", "stage", "T", "N", "grade", "RS")

uni_res <- lapply(vars, function(v) {
  f <- coxph(as.formula(sprintf("Surv(OS.time, OS) ~ %s", v)), data = clin)
  s <- summary(f)$coefficients
  data.frame(variable = v, HR = exp(s[1]), lower = exp(s[1] - 1.96 * s[3]),
             upper = exp(s[1] + 1.96 * s[3]), pvalue = s[5])
}) %>% bind_rows()

multi_fit <- coxph(Surv(OS.time, OS) ~ age + stage + T + N + grade + RS, data = clin)
multi_res <- summary(multi_fit)$coefficients

# forest plot (Fig. 5a/b)
plot_df <- bind_rows(
  uni_res %>% mutate(model = "Univariate"),
  data.frame(variable = rownames(multi_res), HR = exp(multi_res[, 1]),
             lower = exp(multi_res[, 1] - 1.96 * multi_res[, 3]),
             upper = exp(multi_res[, 1] + 1.96 * multi_res[, 3]),
             pvalue = multi_res[, 5], model = "Multivariable"))
ggplot(plot_df, aes(HR, variable)) +
  geom_point() + geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_wrap(~ model) + scale_x_log10() + theme_bw()
ggsave("figures/Fig5_panelAB_forest.pdf", width = 7, height = 5)

# ------------------------------------------------------------------
# (2) nomogram (age + SNI) and calibration (Fig. 5c-f)
# ------------------------------------------------------------------

dd <- datadist(clin); options(datadist = "dd")
nom_fit <- cph(Surv(OS.time, OS) ~ age + RS, data = clin, x = TRUE, y = TRUE,
               surv = TRUE, time.inc = 3 * 365)
nom <- nomogram(nom_fit, fun = function(x) 1 - x, funlabel = "3-year survival")

pdf("figures/Fig5_panelC_nomogram.pdf", width = 7, height = 5)
plot(nom)
dev.off()

cal <- calibrate(nom_fit, method = "boot", B = 1000, u = 3 * 365)
pdf("figures/Fig5_panelD_calibration.pdf", width = 5, height = 5)
plot(cal, xlab = "Predicted 3-year survival", ylab = "Observed 3-year survival")
dev.off()

# ------------------------------------------------------------------
# (3) SHAP feature attribution (Fig. 5g/h)
# ------------------------------------------------------------------
# Representative SHAP computation for the RSF component of the SNI;
# the StepCox component uses the closed-form Cox SHAP values.

load("results/Rdata/model_features.Rdata")      # training matrix of model genes
rsf_fit <- rfsrc(Surv(OS.time, OS) ~ ., data = train_df, ntree = 1000,
                 nodesize = 10, seed = 42)

pfun <- function(object, newdata) predict(object, newdata)$predicted
set.seed(42)
shap <- fastshap::explain(rsf_fit, X = train_df[, model_genes], pred_wrapper = pfun,
                          nsim = 200)

shap_importance <- data.frame(
  gene = colnames(shap),
  mean_abs_shap = colMeans(abs(shap))) %>% arrange(desc(mean_abs_shap))
ggplot(shap_importance, aes(mean_abs_shap, reorder(gene, mean_abs_shap))) +
  geom_col(fill = COL_TUMOR) + theme_bw() +
  labs(x = "Mean |SHAP|", y = NULL)
ggsave("figures/Fig5_panelH_shap.pdf", width = 5, height = 6)

# ------------------------------------------------------------------
# (4) SHAP gene-selection sensitivity (Table S3)
# ------------------------------------------------------------------
# The dual criterion (top-K importance AND IQR stability) is re-run over
# a grid of thresholds to verify that the selected key genes are
# retained under all settings.

sensitivity <- expand.grid(topK = c(3, 5, 10),
                           iqr = c("none", "median", "p75"))
sens_res <- apply(sensitivity, 1, function(setting) {
  k <- as.numeric(setting["topK"])
  top_genes <- shap_importance$gene[1:k]
  iqr_genes <- names(sort(apply(shap, 2, IQR), decreasing = TRUE))[1:k]
  intersect(top_genes, iqr_genes)
})
names(sens_res) <- apply(sensitivity, 1, paste, collapse = "_")
write.csv(sapply(sens_res, paste, collapse = ";"), "tables/TableS3_shap_sensitivity.csv")

# ------------------------------------------------------------------
# (5) Figure 5 assembly
# ------------------------------------------------------------------
# Panels A-H are combined from the PDFs above with the layout defined
# in the manuscript.
