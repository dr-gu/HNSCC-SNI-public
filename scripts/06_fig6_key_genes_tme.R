# ============================================================
# 06_fig6_key_genes_tme.R
# Figure 6 -- key genes (HSPH1 / ST6GALNAC1) and tumour microenvironment
#
# Pipeline:
#   (1) expression validation of the key genes across independent cohorts
#   (2) immune-cell deconvolution of TCGA-HNSC (CIBERSORT)
#   (3) correlation of key-gene expression with immune scores and
#       tumour purity
#   (4) Figure 6 assembly
#
# Requirements: see ENVIRONMENT.md
# ============================================================

rm(list = ls()); gc()

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(corrplot)

COL_TUMOR  <- "#E64B35"
COL_NORMAL <- "#4DBBD5"
KEY_GENES  <- c("HSPH1", "ST6GALNAC1")

# ------------------------------------------------------------------
# (1) expression validation across cohorts (Fig. 6a)
# ------------------------------------------------------------------
# Normalised, batch-corrected matrices for the validation cohorts were
# prepared as described in Methods; tumour/normal labels are in pdata.

cohort_files <- c(
  TCGA = "results/Rdata/expr_TCGA.Rdata",
  `E-MTAB-8588` = "results/Rdata/expr_E_MTAB_8588.Rdata",
  GSE42743 = "results/Rdata/expr_GSE42743.Rdata")

expr_list <- list()
for (nm in names(cohort_files)) {
  load(cohort_files[[nm]])
  expr_list[[nm]] <- expr_corrected[KEY_GENES, , drop = FALSE]
  expr_list[[nm]] <- data.frame(cohort = nm, gene = rownames(expr_corrected[KEY_GENES, ]),
                                t(expr_corrected[KEY_GENES, ]), check.names = FALSE)
}

plot_df <- bind_rows(expr_list) %>%
  pivot_longer(all_of(colnames(expr_list[[1]])[3]), names_to = "sample", values_to = "expr") %>%
  left_join(meta_all, by = c("cohort", "sample"))

ggboxplot(plot_df, x = "group", y = "expr", fill = "group",
          palette = c(COL_NORMAL, COL_TUMOR), facet.by = c("gene", "cohort"),
          scales = "free_y", nrow = 1) +
  stat_compare_means(comparisons = list(c("Normal", "Tumor")),
                     method = "wilcox.test", label = "p.signif", bracket.size = 0) +
  labs(x = NULL, y = "Normalised expression")
ggsave("figures/Fig6_panelA_expression.pdf", width = 12, height = 5)

# ------------------------------------------------------------------
# (2) immune-cell deconvolution (CIBERSORT) (Fig. 6b)
# ------------------------------------------------------------------
# The 22-immune-cell LM22 signature is scored with CIBERSORT
# (third-party implementation under tools/).

source("tools/cibersort.R")
load("data/TCGA_HNSC_expression.Rdata")     # genes x samples TPM matrix

cib <- CIBERSORT(sig_matrix = "data/LM22.txt",
                 mixture_file = "tcga_mixture.txt",
                 perm = 100, QN = TRUE)
# output: sample x 22 cell-fraction matrix (cib)

expr_key <- as.data.frame(t(expr_tpm[KEY_GENES, rownames(cib)]))
colnames(expr_key) <- KEY_GENES

cor_df <- do.call(rbind, lapply(KEY_GENES, function(g) {
  data.frame(gene = g, cell = colnames(cib),
             rho = cor(expr_key[[g]], cib, method = "spearman")[1, ],
             p = apply(cib, 2, function(x) cor.test(expr_key[[g]], x,
                                                    method = "spearman")$p.value))
}))
cor_df$fdr <- p.adjust(cor_df$p, method = "BH")
sig_cor <- cor_df %>% filter(fdr < 0.05) %>%
  mutate(label = ifelse(rho > 0, "Positive", "Negative"))

ggplot(sig_cor, aes(rho, reorder(cell, rho), fill = label)) +
  geom_col() + facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c(Positive = COL_TUMOR, Negative = COL_NORMAL)) +
  theme_bw() + labs(x = "Spearman rho", y = NULL)
ggsave("figures/Fig6_panelB_deconvolution.pdf", width = 8, height = 6)

# ------------------------------------------------------------------
# (3) correlation with ESTIMATE scores (Fig. 6c)
# ------------------------------------------------------------------

load("results/estimate_scores.Rdata")          # immune / stromal / purity
est_df <- est_df[rownames(expr_key), ]
cor_est <- data.frame(gene = rep(KEY_GENES, each = ncol(est_df)),
                      score = rep(colnames(est_df), length(KEY_GENES)),
                      rho = c(cor(expr_key[, 1], est_df, method = "spearman")[1, ],
                              cor(expr_key[, 2], est_df, method = "spearman")[1, ]))
ggplot(cor_est, aes(gene, score, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = COL_NORMAL, mid = "white", high = COL_TUMOR) +
  geom_text(aes(label = round(rho, 2))) + theme_minimal() +
  labs(x = NULL, y = NULL, fill = "rho")
ggsave("figures/Fig6_panelC_estimate_cor.pdf", width = 5, height = 3.5)

# ------------------------------------------------------------------
# (4) Figure 6 assembly
# ------------------------------------------------------------------
# Panels A-C are combined from the PDFs above with the layout defined
# in the manuscript.
