library("here")
library("SummarizedExperiment")
library("readxl")
library("jaffelab")
library("edgeR")
library("sessioninfo")

## For reproducing the jitter output
set.seed(20230907)

## Create output directories
dir_plots <- here("plots", "06_check_gene_exprs")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir_rdata <- here("processed-data", "06_check_gene_exprs")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## published results
res <-
    read_excel(here(
        "raw-data",
        "LBP",
        "DONE_LEL2023_LBP_FLAGSHIP_SCIENCE_ST1.xlsx"
    ),
        sheet = 2)
res <- res[order(res$P.Value), ]
colnames(res)[c(1, 8)] <- c("ensembl_id", "de_status")

## Load SPEAQeasy gene level data
rse_gene <- readRDS(here("processed-data", "02_SPEAQeasy", "rse_gene_living_brain_reanalysis_n516.Rds"))
rowData(rse_gene)$ensembl_id <- rowData(rse_gene)$ensemblID

## phenotype data
rse_gene$COI <- factor(ifelse(rse_gene$isPostMortem, "PM", "LIV"))
rse_gene$postmortem <- as.numeric(rse_gene$COI) - 1

## add cell comp PCs
cellProps <- read.csv(here("processed-data", "04_deconvolution", "LBP_burkeDecon.csv"), row.names = 1)
cellProps <- cellProps[colnames(rse_gene), ]
cellPca <- prcomp(cellProps)
getPcaVars(cellPca)[1:2]
# [1] 79.0 11.1
rse_gene$cellPC <- cellPca$x[, 1]

write.csv(as.data.frame(colData(rse_gene)),
    file = file.path(dir_rdata, "LBP_merged_seqmetrics.csv"),
    row.names = FALSE
)

## normalize
dge <- DGEList(
    counts = assays(rse_gene)$counts,
    genes = rowData(rse_gene)
)
dge <- calcNormFactors(dge)

### PCA
log_rpkm <- rpkm(dge, gene.length = rowData(rse_gene)$score, log = TRUE)
pca <- prcomp(t(log_rpkm))

##############
## modeling ##
##############

# use their genes for FDR calc
eIndex <- which(rowData(rse_gene)$ensembl_id %in% res$ensembl_id)

## For simplicity later on, I'll use this variable name
rse_gene$samplingAge <- rse_gene$ageDeath_num.biospecimen

## Fix race since it has 206 NAs
rse_gene$race <- as.character(rse_gene$race)
rse_gene$race[is.na(rse_gene$race)] <- "Unknown"

## model
mod <- model.matrix(
    ~ COI + race + sex + diagnosis +
        RIN + cellPC +
        mitoRate +
        totalAssignedGene +
        rRNA_rate +
        overallMapRate,
    data = colData(rse_gene)
)
qsvs_cell <- read.csv(here("processed-data", "05_check_degradation", "qSVs_cell.csv"), row.names = 1)
mod_qsva <- cbind(mod, qsvs_cell[rownames(mod), ])

mod_pc3 <- cbind(mod, PC3 = pca$x[, 3])

## voom
vGene <- voom(dge, mod, plot = FALSE)
vGene_qsva <- voom(dge, mod_qsva, plot = FALSE)
vGene_pc3 <- voom(dge, mod_pc3, plot = FALSE)

eBGene <- eBayes(lmFit(vGene))
eBGene_qsva <- eBayes(lmFit(vGene_qsva))
eBGene_pc3 <- eBayes(lmFit(vGene_pc3))

## extract
outGene <- topTable(
    eBGene,
    coef = 2,
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene$isExprs <- outGene$ensembl_id %in% res$ensembl_id
table(outGene$isExprs)
# FALSE  TRUE
# 36882 21155
dim(outGene)
# [1] 58037    18
dim(res)
# [1] 21635     8
outGene$adj.P.Val <- NA
outGene$adj.P.Val[outGene$isExprs] <- p.adjust(outGene$P.Value[outGene$isExprs], "fdr")

outGene_qsva <- topTable(
    eBGene_qsva,
    coef = 2,
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_qsva$isExprs <- outGene_qsva$ensembl_id %in% res$ensembl_id
outGene_qsva$adj.P.Val <- NA
outGene_qsva$adj.P.Val[outGene_qsva$isExprs] <- p.adjust(outGene_qsva$P.Value[outGene_qsva$isExprs], "fdr")

outGene_pc3 <- topTable(
    eBGene_pc3,
    coef = 2,
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_pc3$isExprs <- outGene_pc3$ensembl_id %in% res$ensembl_id
outGene_pc3$adj.P.Val <- NA
outGene_pc3$adj.P.Val[outGene_pc3$isExprs] <- p.adjust(outGene_pc3$P.Value[outGene_pc3$isExprs], "fdr")

table(outGene$adj.P.Val < 0.05)
# FALSE  TRUE
#  5602 15553
table(outGene_qsva$adj.P.Val < 0.05)
# FALSE  TRUE
# 17423  3732
table(outGene_pc3$adj.P.Val < 0.05)
# FALSE  TRUE
#  5736 15419

outGene_order <- outGene[order(outGene$P.Value), ]
outGene_order$mean_dir <- cumsum(sign(outGene_order$t)) / 1:nrow(outGene_order)

pdf(file.path(dir_plots, "outGene_order_mean_dir.pdf"))
plot(outGene_order$mean_dir)
dev.off()

pdf(file.path(dir_plots, "published_vs_redone_tstat_modelCompare.pdf"))
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 2,
    cex.lab = 2
)
res_match <- res[match(outGene$ensembl_id, res$ensembl_id), ]
plot(
    res_match$t,
    outGene$t,
    xlab = "Pre-print DE t-stat",
    ylab = "Re-calculated DE t-stat",
    pch = 21,
    bg = "grey"
)
abline(0, 1, lty = 2, col = "red")
dev.off()

pdf(file.path(dir_plots, "qsva_attentuation.pdf"))
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 2,
    cex.lab = 2
)
plot(
    outGene$t[outGene$isExprs] ~ outGene_qsva$t[outGene_qsva$isExprs],
    pch = 21,
    bg = "grey",
    xlab = "PM vs LIV t-stat (with qSVA)",
    ylab = "PM vs LIV t-stat (without qSVA)"
)
abline(0,
    1,
    col = "blue",
    lwd = 3,
    lty = 2
)

plot(
    outGene$logFC[outGene$isExprs] ~ outGene_qsva$logFC[outGene_qsva$isExprs],
    pch = 21,
    bg = "grey",
    xlab = "PM vs LIV log2FC (with qSVA)",
    ylab = "PM vs LIV log2FC (without qSVA)"
)
abline(0,
    1,
    col = "blue",
    lwd = 3,
    lty = 2
)
dev.off()

###########
## age in liv
rse_liv <- rse_gene[, rse_gene$COI == "LIV"]
table(rse_liv$diagnosis)
rse_liv$diagnosis <- droplevels(rse_liv$diagnosis)
table(rse_liv$diagnosis)

## model
mod_liv <- model.matrix(
    ~ samplingAge + race + sex + diagnosis +
        RIN + cellPC +
        mitoRate +
        totalAssignedGene +
        rRNA_rate +
        overallMapRate,
    data = colData(rse_liv)
)
qsvs_liv_cell <- read.csv(here("processed-data", "05_check_degradation", "qSVs_cell_liv.csv"), row.names = 1)
mod_liv_qsva <- cbind(mod_liv, qsvs_liv_cell[rownames(mod_liv), ])
dge_liv <- dge[, rownames(mod_liv)]

## voom
vGene_liv <- voom(dge_liv, mod_liv, plot = FALSE)
vGene_liv_qsva <- voom(dge_liv, mod_liv_qsva, plot = FALSE)

eBGene_liv <- eBayes(lmFit(vGene_liv))
eBGene_liv_qsva <- eBayes(lmFit(vGene_liv_qsva))

## extract
outGene_liv <- topTable(
    eBGene_liv,
    coef = "samplingAge",
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_liv$isExprs <- outGene_liv$ensembl_id %in% res$ensembl_id
outGene_liv$adj.P.Val <- NA
outGene_liv$adj.P.Val[outGene_liv$isExprs] <- p.adjust(outGene_liv$P.Value[outGene_liv$isExprs], "fdr")

outGene_qsva_liv <- topTable(
    eBGene_liv_qsva,
    coef = "samplingAge",
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_qsva_liv$isExprs <- outGene_qsva_liv$ensembl_id %in% res$ensembl_id
outGene_qsva_liv$adj.P.Val <- NA
outGene_qsva_liv$adj.P.Val[outGene_qsva_liv$isExprs] <- p.adjust(outGene_qsva_liv$P.Value[outGene_qsva_liv$isExprs], "fdr")

table(outGene_liv$adj.P.Val < 0.05)
# FALSE  TRUE
# 20609   546
table(outGene_qsva_liv$adj.P.Val < 0.05)
# FALSE  TRUE
# 19330  1825

###########
## age in pm
rse_pm <- rse_gene[, rse_gene$COI == "PM"]
table(rse_pm$diagnosis)
rse_pm$diagnosis <- droplevels(rse_pm$diagnosis)
table(rse_pm$diagnosis)

## model
mod_pm <- model.matrix(
    ~ samplingAge + race + sex + diagnosis +
        RIN + cellPC +
        mitoRate +
        totalAssignedGene +
        rRNA_rate +
        overallMapRate,
    data = colData(rse_pm)
)

qsvs_pm_cell <- read.csv(here("processed-data", "05_check_degradation", "qSVs_cell_pm.csv"), row.names = 1)
mod_pm_qsva <- cbind(mod_pm, qsvs_pm_cell[rownames(mod_pm), ])
dge_pm <- dge[, rownames(mod_pm)]

## voom
vGene_pm <- voom(dge_pm, mod_pm, plot = FALSE)
vGene_pm_qsva <- voom(dge_pm, mod_pm_qsva, plot = FALSE)

eBGene_pm <- eBayes(lmFit(vGene_pm))
eBGene_pm_qsva <- eBayes(lmFit(vGene_pm_qsva))

## extract
outGene_pm <- topTable(
    eBGene_pm,
    coef = "samplingAge",
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_pm$isExprs <- outGene_pm$ensembl_id %in% res$ensembl_id
outGene_pm$adj.P.Val <- NA
outGene_pm$adj.P.Val[outGene_pm$isExprs] <- p.adjust(outGene_pm$P.Value[outGene_pm$isExprs], "fdr")

outGene_qsva_pm <- topTable(
    eBGene_pm_qsva,
    coef = "samplingAge",
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_qsva_pm$isExprs <- outGene_qsva_pm$ensembl_id %in% res$ensembl_id
outGene_qsva_pm$adj.P.Val <- NA
outGene_qsva_pm$adj.P.Val[outGene_qsva_pm$isExprs] <- p.adjust(outGene_qsva_pm$P.Value[outGene_qsva_pm$isExprs], "fdr")

table(outGene_pm$adj.P.Val < 0.05)
# FALSE  TRUE
# 16920  4235
table(outGene_qsva_pm$adj.P.Val < 0.05)
# FALSE  TRUE
# 19788  1367

pdf(file.path(dir_plots, "age_pm_logFC_vs_age_liv_logFC.pdf"))
plot(
    outGene_pm$logFC[outGene_pm$isExprs],
    outGene_liv$logFC[outGene_pm$isExprs],
    xlab = "Age logFC (PM, no qSVA)",
    ylab = "Age logFC (LIV, no qSVA)",
    pch = 21,
    bg = "grey"
)
dev.off()

cor.test(outGene_pm$logFC[outGene_pm$isExprs], outGene_liv$logFC[outGene_pm$isExprs])
# 	Pearson's product-moment correlation
#
# data:  outGene_pm$logFC[outGene_pm$isExprs] and outGene_liv$logFC[outGene_pm$isExprs]
# t = 36.316, df = 21153, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2295303 0.2548999
# sample estimates:
#       cor
# 0.2422565

pdf(file.path(dir_plots, "age_qsva_pm_logFC_vs_age_qsva_liv_logFC.pdf"))
plot(
    outGene_qsva_pm$logFC[outGene_pm$isExprs],
    outGene_qsva_liv$logFC[outGene_pm$isExprs],
    xlab = "Age logFC (PM, qSVA)",
    ylab = "Age logFC (LIV, qSVA)",
    pch = 21,
    bg = "grey"
)
dev.off()
cor.test(outGene_qsva_pm$logFC[outGene_pm$isExprs], outGene_qsva_liv$logFC[outGene_pm$isExprs])
# 	Pearson's product-moment correlation
#
# data:  outGene_qsva_pm$logFC[outGene_pm$isExprs] and outGene_qsva_liv$logFC[outGene_pm$isExprs]
# t = 43.731, df = 21153, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2755385 0.3002554
# sample estimates:
#       cor
# 0.2879449

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 (2023-06-16)
#  os       macOS Ventura 13.5
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2023-09-08
#  rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
#  pandoc   3.1.5 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
#  annotate               1.78.0    2023-04-25 [1] Bioconductor
#  AnnotationDbi          1.62.2    2023-07-02 [1] Bioconductor
#  Biobase              * 2.60.0    2023-04-25 [1] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [1] Bioconductor
#  BiocParallel           1.34.2    2023-05-28 [1] Bioconductor
#  Biostrings             2.68.1    2023-05-16 [1] Bioconductor
#  bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
#  bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
#  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  blob                   1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
#  brio                   1.1.3     2021-11-30 [1] CRAN (R 4.3.0)
#  cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  callr                  3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.1)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  data.table             1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
#  DBI                    1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
#  DelayedArray           0.26.7    2023-07-30 [1] Bioconductor
#  devtools             * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
#  digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
#  dplyr                  1.1.3     2023-09-03 [1] CRAN (R 4.3.0)
#  edgeR                * 3.42.4    2023-05-31 [1] Bioconductor
#  ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
#  fansi                  1.0.4     2023-01-22 [1] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
#  fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
#  gargle                 1.5.2     2023-07-20 [1] CRAN (R 4.3.0)
#  genefilter             1.82.1    2023-05-02 [1] Bioconductor
#  generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.2    2023-08-27 [1] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-05-06 [1] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [1] Bioconductor
#  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
#  googledrive            2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
#  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools              0.5.6     2023-08-10 [1] CRAN (R 4.3.0)
#  htmlwidgets            1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
#  httpuv                 1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
#  httr                   1.4.7     2023-08-15 [1] CRAN (R 4.3.0)
#  IRanges              * 2.34.1    2023-07-02 [1] Bioconductor
#  jaffelab             * 0.99.32   2023-05-06 [1] Github (LieberInstitute/jaffelab@7b7afe3)
#  KEGGREST               1.40.0    2023-04-25 [1] Bioconductor
#  later                  1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
#  lattice                0.21-8    2023-04-05 [1] CRAN (R 4.3.1)
#  lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
#  limma                * 3.56.2    2023-06-04 [1] Bioconductor
#  locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
#  lubridate              1.9.2     2023-02-10 [1] CRAN (R 4.3.0)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
#  MASS                   7.3-60    2023-05-04 [1] CRAN (R 4.3.1)
#  Matrix                 1.6-1     2023-08-14 [1] CRAN (R 4.3.0)
#  MatrixGenerics       * 1.12.3    2023-07-30 [1] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
#  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
#  mgcv                   1.9-0     2023-07-11 [1] CRAN (R 4.3.0)
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
#  nlme                   3.1-163   2023-08-09 [1] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild               1.4.2     2023-06-26 [1] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload                1.3.2.1   2023-07-08 [1] CRAN (R 4.3.0)
#  png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
#  prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
#  processx               3.8.2     2023-06-30 [1] CRAN (R 4.3.0)
#  profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
#  promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
#  prompt                 1.0.2     2023-08-31 [1] CRAN (R 4.3.0)
#  ps                     1.7.5     2023-04-18 [1] CRAN (R 4.3.0)
#  purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
#  RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
#  RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
#  readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
#  remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
#  rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
#  RSQLite                2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
#  rsthemes               0.4.0     2023-05-06 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
#  S4Arrays               1.0.5     2023-07-30 [1] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [1] Bioconductor
#  segmented              1.6-4     2023-04-13 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
#  shiny                  1.7.5     2023-08-12 [1] CRAN (R 4.3.0)
#  stringi                1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
#  stringr                1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
#  SummarizedExperiment * 1.30.2    2023-06-11 [1] Bioconductor
#  suncalc                0.5.1     2022-09-29 [1] CRAN (R 4.3.0)
#  survival               3.5-7     2023-08-14 [1] CRAN (R 4.3.0)
#  sva                    3.48.0    2023-04-25 [1] Bioconductor
#  testthat             * 3.1.10    2023-07-06 [1] CRAN (R 4.3.0)
#  tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
#  timechange             0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
#  urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
#  usethis              * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
#  vctrs                  0.6.3     2023-06-14 [1] CRAN (R 4.3.0)
#  XML                    3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-04-25 [1] Bioconductor
#  zlibbioc               1.46.0    2023-04-25 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
