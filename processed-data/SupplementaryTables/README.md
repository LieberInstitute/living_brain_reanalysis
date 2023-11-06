Supplementary Tables
====================

## Table S1

Table with the p-values for the correlation between the top 10 gene-level principal components (PCs) and several _SPEAQeasy_ metrics as well as the LIV vs PM indicator variable. The -log10 p-values are displayed in **Figure 1A**. Created [here](https://github.com/LieberInstitute/living_brain_reanalysis/blob/dbe7d06ff7bc148eb29a9fac11f5839d5a1a31ce/code/03_check_QC/01_qc_plots.R#L98).

## Table S2

Merged sequencing metrics produced by _SPEAQeasy_, first principal component from the estimated cell type proportions (_cellPC_), and other sample demographic information. See <https://research.libd.org/SPEAQeasy/outputs.html#quality-metrics> for details about the _SPEAQeasy_ alignment metrics. Data from this table is used for **Figure 1B**. Created [here](https://github.com/LieberInstitute/living_brain_reanalysis/blob/dbe7d06ff7bc148eb29a9fac11f5839d5a1a31ce/code/06_check_gene_exprs/01_check_gene_exprs.R#L45).

## Table S3

Per sample, estimated cell type proportions via deconvolution. Data from this table is used for **Figure 1C**. Created [here](https://github.com/LieberInstitute/living_brain_reanalysis/blob/c856e4712d469630cd5618d136f17477dd0275ae/code/04_deconvolution/01_deconvolution_burkeModel.R#L49).


## Table S4

Limma-voom gene-level differential expression results between LIV and PM samples either (a) without adjusting for qSVs (quality Surrogate Variables), (b) adjusting for qSVs, or (c) reproducing the original results by Liharska et al. Data from this table is used for **Figure 1D** and **Figure 1F**.Created [here](https://github.com/LieberInstitute/living_brain_reanalysis/blob/dbe7d06ff7bc148eb29a9fac11f5839d5a1a31ce/code/06_check_gene_exprs/01_check_gene_exprs.R#L233).

## Table S5

Limma transcript-level differential expression results between (a) LIV and PM samples or (b) by RNA degradation (Jaffe et al., _PNAS_, 2017). Data from this table is used for **Figure 1E**.Created [here](https://github.com/LieberInstitute/living_brain_reanalysis/blob/3a11b296b70b6b8a81911d9e5800cc5f9f6b889f/code/05_check_degradation/02_check_deg_matrix.R#L568).
