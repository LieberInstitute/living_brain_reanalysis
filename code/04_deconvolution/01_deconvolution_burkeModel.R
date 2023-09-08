library("readxl")
library("jaffelab")
library("RColorBrewer")
library("here")
library("SummarizedExperiment")
library("edgeR")
library("minfi")
library("sessioninfo")

## For reproducing the jitter output
set.seed(20230907)

## Create output directories
dir_plots <- here("plots", "04_deconvolution")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir_rdata <- here("processed-data", "04_deconvolution")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## read in decon
decon_df <- read_excel(here("raw-data", "Burke", "decon_model.xlsx"), sheet = 1, skip = 3)
colnames(decon_df)[3] <- "CellType"

## make matrix
coefEsts <- as.matrix(decon_df[, 4:13])
rownames(coefEsts) <- ss(decon_df$Gene, "\\.")

## Load SPEAQeasy gene level data
rse_gene <- readRDS(here("processed-data", "02_SPEAQeasy", "rse_gene_living_brain_reanalysis_n516.Rds"))
rowData(rse_gene)$ensembl_id <- rowData(rse_gene)$ensemblID

## phenotype data
rse_gene$COI <- factor(ifelse(rse_gene$isPostMortem, "PM", "LIV"))
rse_gene$postmortem <- as.numeric(rse_gene$COI) - 1

## get RPKM
dge <- DGEList(
    counts = assays(rse_gene)$counts,
    genes = rowData(rse_gene)
)
dge <- calcNormFactors(dge)
yExprs <- rpkm(dge, gene.length = rowData(rse_gene)$Length, log = TRUE)

rownames(yExprs) <- rowData(rse_gene)$ensembl_id
yExprs_scaled <- scale(yExprs[rownames(coefEsts), ]) # filter to cell type genes, and scale

## do deconvolution
propEsts <- minfi:::projectCellType(yExprs_scaled, coefEsts)
propEsts_scaled <- prop.table(propEsts, 1)
write.csv(propEsts_scaled, file = file.path(dir_rdata, "LBP_burkeDecon.csv"))

##  each variable vs COI
pdf(file.path(dir_plots, "cellTypes_vs_COI.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 2,
    cex.lab = 2
)
palette(brewer.pal(4, "Dark2"))
for (i in 1:ncol(propEsts_scaled)) {
    boxplot(
        propEsts_scaled[, i] ~ rse_gene$COI,
        ylab = colnames(propEsts_scaled)[i],
        xlab = "",
        outline = FALSE,
        ylim = range(propEsts_scaled[, i])
    )
    points(
        propEsts_scaled[, i] ~
            jitter(as.numeric(rse_gene$COI), amount = 0.15),
        pch = 21,
        bg = rse_gene$COI
    )
    tt <- summary(lm(propEsts_scaled[, i] ~ rse_gene$COI))$coef[2, ]
    lc <- ifelse(tt[1] > 0, "topleft", "topright")
    legend(lc, paste0("p=", signif(tt[4], 3)), cex = 1.5)
}
dev.off()

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
#  date     2023-09-07
#  rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
#  pandoc   3.1.5 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
#  annotate               1.78.0    2023-04-25 [1] Bioconductor
#  AnnotationDbi          1.62.2    2023-07-02 [1] Bioconductor
#  askpass                1.2.0     2023-09-03 [1] CRAN (R 4.3.0)
#  base64                 2.0.1     2022-08-19 [1] CRAN (R 4.3.0)
#  beanplot               1.3.1     2022-04-09 [1] CRAN (R 4.3.0)
#  Biobase              * 2.60.0    2023-04-25 [1] Bioconductor
#  BiocFileCache          2.8.0     2023-04-25 [1] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [1] Bioconductor
#  BiocIO                 1.10.0    2023-04-25 [1] Bioconductor
#  BiocParallel           1.34.2    2023-05-28 [1] Bioconductor
#  biomaRt                2.56.1    2023-06-11 [1] Bioconductor
#  Biostrings           * 2.68.1    2023-05-16 [1] Bioconductor
#  bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
#  bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
#  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  blob                   1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
#  brio                   1.1.3     2021-11-30 [1] CRAN (R 4.3.0)
#  bumphunter           * 1.42.0    2023-04-25 [1] Bioconductor
#  cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  callr                  3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.1)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  curl                   5.0.2     2023-08-14 [1] CRAN (R 4.3.1)
#  data.table             1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
#  DBI                    1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
#  dbplyr                 2.3.3     2023-07-07 [1] CRAN (R 4.3.0)
#  DelayedArray           0.26.7    2023-07-30 [1] Bioconductor
#  DelayedMatrixStats     1.22.5    2023-08-13 [1] Bioconductor
#  devtools             * 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
#  digest                 0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
#  doRNG                  1.8.6     2023-01-16 [1] CRAN (R 4.3.0)
#  dplyr                  1.1.3     2023-09-03 [1] CRAN (R 4.3.0)
#  edgeR                * 3.42.4    2023-05-31 [1] Bioconductor
#  ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
#  fansi                  1.0.4     2023-01-22 [1] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
#  filelock               1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
#  foreach              * 1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
#  fs                     1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
#  gargle                 1.5.2     2023-07-20 [1] CRAN (R 4.3.0)
#  genefilter             1.82.1    2023-05-02 [1] Bioconductor
#  generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.2    2023-08-27 [1] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-05-06 [1] Bioconductor
#  GenomicAlignments      1.36.0    2023-04-25 [1] Bioconductor
#  GenomicFeatures        1.52.2    2023-08-27 [1] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [1] Bioconductor
#  GEOquery               2.68.0    2023-04-25 [1] Bioconductor
#  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
#  googledrive            2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
#  HDF5Array              1.28.1    2023-05-01 [1] Bioconductor
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
#  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools              0.5.6     2023-08-10 [1] CRAN (R 4.3.0)
#  htmlwidgets            1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
#  httpuv                 1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
#  httr                   1.4.7     2023-08-15 [1] CRAN (R 4.3.0)
#  illuminaio             0.42.0    2023-05-08 [1] Bioconductor
#  IRanges              * 2.34.1    2023-07-02 [1] Bioconductor
#  iterators            * 1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
#  jaffelab             * 0.99.32   2023-05-06 [1] Github (LieberInstitute/jaffelab@7b7afe3)
#  KEGGREST               1.40.0    2023-04-25 [1] Bioconductor
#  later                  1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
#  lattice                0.21-8    2023-04-05 [1] CRAN (R 4.3.1)
#  lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
#  limma                * 3.56.2    2023-06-04 [1] Bioconductor
#  locfit               * 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
#  lubridate              1.9.2     2023-02-10 [1] CRAN (R 4.3.0)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
#  MASS                   7.3-60    2023-05-04 [1] CRAN (R 4.3.1)
#  Matrix                 1.6-1     2023-08-14 [1] CRAN (R 4.3.0)
#  MatrixGenerics       * 1.12.3    2023-07-30 [1] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
#  mclust                 6.0.0     2022-10-31 [1] CRAN (R 4.3.0)
#  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  minfi                * 1.46.0    2023-05-08 [1] Bioconductor
#  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
#  multtest               2.56.0    2023-05-08 [1] Bioconductor
#  nlme                   3.1-163   2023-08-09 [1] CRAN (R 4.3.0)
#  nor1mix                1.3-0     2019-06-13 [1] CRAN (R 4.3.0)
#  openssl                2.1.0     2023-07-15 [1] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild               1.4.2     2023-06-26 [1] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload                1.3.2.1   2023-07-08 [1] CRAN (R 4.3.0)
#  plyr                   1.8.8     2022-11-11 [1] CRAN (R 4.3.0)
#  png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
#  preprocessCore         1.62.1    2023-05-08 [1] Bioconductor
#  prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
#  processx               3.8.2     2023-06-30 [1] CRAN (R 4.3.0)
#  profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
#  progress               1.2.2     2019-05-16 [1] CRAN (R 4.3.0)
#  promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
#  prompt                 1.0.2     2023-08-31 [1] CRAN (R 4.3.0)
#  ps                     1.7.5     2023-04-18 [1] CRAN (R 4.3.0)
#  purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
#  quadprog               1.5-8     2019-11-20 [1] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
#  rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
#  RColorBrewer         * 1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
#  RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
#  readr                  2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
#  readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
#  remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
#  reshape                0.8.9     2022-04-12 [1] CRAN (R 4.3.0)
#  restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.3.0)
#  rhdf5                  2.44.0    2023-04-25 [1] Bioconductor
#  rhdf5filters           1.12.1    2023-04-30 [1] Bioconductor
#  Rhdf5lib               1.22.0    2023-04-25 [1] Bioconductor
#  rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
#  rngtools               1.5.2     2021-09-20 [1] CRAN (R 4.3.0)
#  rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
#  Rsamtools              2.16.0    2023-04-25 [1] Bioconductor
#  RSQLite                2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
#  rsthemes               0.4.0     2023-05-06 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
#  rtracklayer            1.60.1    2023-08-20 [1] Bioconductor
#  S4Arrays               1.0.5     2023-07-30 [1] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [1] Bioconductor
#  scrime                 1.3.5     2018-12-01 [1] CRAN (R 4.3.0)
#  segmented              1.6-4     2023-04-13 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
#  shiny                  1.7.5     2023-08-12 [1] CRAN (R 4.3.0)
#  siggenes               1.74.0    2023-05-08 [1] Bioconductor
#  sparseMatrixStats      1.12.2    2023-07-02 [1] Bioconductor
#  stringi                1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
#  stringr                1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
#  SummarizedExperiment * 1.30.2    2023-06-11 [1] Bioconductor
#  suncalc                0.5.1     2022-09-29 [1] CRAN (R 4.3.0)
#  survival               3.5-7     2023-08-14 [1] CRAN (R 4.3.0)
#  testthat             * 3.1.10    2023-07-06 [1] CRAN (R 4.3.0)
#  tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
#  tidyr                  1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
#  timechange             0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
#  tzdb                   0.4.0     2023-05-12 [1] CRAN (R 4.3.0)
#  urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
#  usethis              * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
#  vctrs                  0.6.3     2023-06-14 [1] CRAN (R 4.3.0)
#  XML                    3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
#  xml2                   1.3.5     2023-07-06 [1] CRAN (R 4.3.0)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
#  XVector              * 0.40.0    2023-04-25 [1] Bioconductor
#  yaml                   2.3.7     2023-01-23 [1] CRAN (R 4.3.0)
#  zlibbioc               1.46.0    2023-04-25 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
