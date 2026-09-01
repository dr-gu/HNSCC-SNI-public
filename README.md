# HNSCC-SNI-public

Analysis code accompanying the HNSC-SNI study: construction and
validation of a Sialylation-Niche Index (SNI) prognostic model for head
and neck squamous cell carcinoma (HNSCC), including molecular subtyping,
single-cell and spatial characterisation of the key genes HSPH1 and
ST6GALNAC1, and an AI-guided drug-screening pipeline.

## Repository contents

| Path | Description |
|---|---|
| `scripts/` | One analysis script per figure, following the manuscript's figure structure |
| `screening/` | Representative excerpt of the multi-stage drug-screening pipeline (GraphBAN / ADMET-AI / docking / GROMACS) |
| `ENVIRONMENT.md` | Machine, R/Python versions and key packages used |
| `LICENSE` | MIT license |

### Scripts

| Script | Figure | Analysis |
|---|---|---|
| `01_fig2_sialo_deg_landscape.R` | Fig. 2 | Differential expression across GEO cohorts, sialylation-gene intersection, univariate Cox screen, GO/KEGG enrichment, mutation landscape |
| `02_fig3_molecular_subtypes.R` | Fig. 3 | Consensus clustering (C1/C2), survival, GSEA, ESTIMATE, immune checkpoints, TIP steps |
| `03_fig4_sni_construction.R` | Fig. 4 (construction) | 101-combination machine-learning benchmark (Mime1), model selection, risk score, KM, time-dependent ROC |
| `04_fig4_sni_validation_robustness.R` | Fig. 4 (validation); Figs. S1–S3 | External-cohort validation, bootstrap CI, pooled validation, clinical benchmark |
| `05_fig5_nomogram_shap.R` | Fig. 5 | Multivariable Cox, nomogram + calibration, SHAP attribution, sensitivity grid |
| `06_fig6_key_genes_tme.R` | Fig. 6 | Key-gene expression validation, CIBERSORT deconvolution, immune-score correlations |
| `07_fig7_scRNA_knockout.R` | Fig. 7 | Seurat single-cell workflow, cell-type distribution, scTenifoldKnk virtual knockout |
| `08_fig8_spatial.R` | Fig. 8 | 10x Visium spatial transcriptomics, niche scoring, colocalisation |
| `09_fig9_drug_sensitivity_md.R` | Fig. 9 | Pharmacogenomic correlations (GDSC2 / CTRP / PRISM / oncoPredict bridge) |
| `09_fig9_md_analysis.py` | Fig. 9 | MD trajectory metrics (RMSD, contacts, minimum distance, hydrogen bonds) |

## Data availability

The scripts read processed matrices and annotation files that are
obtained from public repositories as described in the manuscript
Methods: GEO (GSE29330, GSE30784, GSE3292, GSE42743, GSE6791, GSE7224,
GSE9844, GSE65858, GSE41613, GSE234933, GSE181300), ArrayExpress
(E-MTAB-8588), TCGA (TCGA-HNSC expression, clinical, and mutation data),
GDSC2, CTRP v2 and DepMap PRISM. Intermediate processed files are not
redistributed with this repository.

## Citation

Please cite the published article.
