library("readxl")
library("jaffelab")
library("edgeR")
library("here")
library("SummarizedExperiment")
library("sessioninfo")

## For reproducing the jitter output
set.seed(20230907)

pcHeatmap <-
    function(rse,
    pca,
    plotdir,
    pthresh = 10,
    width = 12,
    height = 14) {
        library("RColorBrewer")
        library("lattice")
        n <- min(10, ncol(pca$x))
        vars <- colnames(colData(rse))[sapply(colData(rse), class) ==
            "numeric"]
        vars <- vars[!grepl("ageDeath_num.individual|_R2|_R1|readLength", vars)]
        ccPcaMetrics <- cor(pca$x[, 1:n], as.data.frame(colData(rse)[
            ,
            vars
        ]), use = "complete.obs")
        tPcaMetrics <-
            ccPcaMetrics / sqrt((1 - ccPcaMetrics^2) / (ncol(rse) - 2))
        pvalPcaMetrics <- 2 * pt(-abs(tPcaMetrics), df = ncol(rse) - 1)
        filepath <- file.path(plotdir, "pcs_vs_metrics_heatmap.pdf")
        filepath <- gsub("//", "/", filepath)
        pdf(filepath, w = width, h = height)
        logPvalMat <- -log10(pvalPcaMetrics)
        logPvalMat[logPvalMat > pthresh] <- pthresh
        theSeq <- seq(0, pthresh, by = 0.1)
        my.col <- colorRampPalette(c("white", "blue"))(length(theSeq))
        print(
            levelplot(
                logPvalMat,
                aspect = "fill",
                at = theSeq,
                pretty = TRUE,
                xlab = "",
                ylab = "",
                scales = list(x = list(
                    rot = 90,
                    cex = 1.2
                ), y = list(cex = 1.2)),
                panel = panel.levelplot.raster,
                col.regions = my.col
            )
        )
        dev.off()
        return(pvalPcaMetrics)
    }

## Create output directory for plots
dir_plots <- here("plots", "03_check_QC")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

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

### do PCA
dge <- DGEList(
    counts = assays(rse_gene)$counts,
    genes = rowData(rse_gene)
)
dge <- calcNormFactors(dge)

### PCA
log_rpkm <-
    rpkm(dge, gene.length = rowData(rse_gene)$Length, log = TRUE)
pca <- prcomp(t(log_rpkm))
pcaVars <- getPcaVars(pca)
pcaVars[1:10]
# [1] 27.600  6.150  4.610  2.870  2.360  1.920  1.780  1.440  1.110  0.919

pcMat <- pcHeatmap(rse_gene, pca, plotdir = dir_plots, pthresh = 20)
dir.create(here("processed-data", "SupplementaryTables"))
write.csv(pcMat, file = here("processed-data", "SupplementaryTables", "TableS1.csv"))

###### example plots
vars <-
    colnames(colData(rse_gene))[sapply(colData(rse_gene), class) == "numeric"]
vars <- vars[!grepl("ageDeath_num.individual|_R2|_R1|readLength|postmortem", vars)]
ccPcaMetrics <-
    cor(pca$x[, 1:10], as.data.frame(colData(rse_gene)[, vars]), use = "complete.obs")
nn <- c("totalAssignedGene")
ccPcaMetrics[1:3, nn]
gIndexes <- splitit(rse_gene$COI)

##  each variable vs COI
pdf(file.path(dir_plots, "Vars_vs_COI.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 2,
    cex.lab = 2
)
palette(brewer.pal(4, "Dark2"))
for (i in seq(along = vars)) {
    boxplot(
        colData(rse_gene)[, vars[i]] ~ rse_gene$COI,
        ylab = vars[i],
        xlab = "",
        outline = FALSE,
        ylim = range(colData(rse_gene)[, vars[i]], na.rm = TRUE)
    )
    points(
        colData(rse_gene)[, vars[i]] ~
            jitter(as.numeric(rse_gene$COI), amount = 0.15),
        pch = 21,
        bg = rse_gene$COI
    )
    tt <-
        summary(lm(colData(rse_gene)[, vars[i]] ~ rse_gene$COI))$coef[2, ]
    lc <- ifelse(tt[1] > 0, "topleft", "topright")
    legend(lc, paste0("p=", signif(tt[4], 3)), cex = 1.5)
}
dev.off()


## PC1
pdf(file.path(dir_plots, "PC1_vs_COI.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
boxplot(
    pca$x[, 1] ~ rse_gene$COI,
    xlab = "",
    ylab = paste0("PC1: ", pcaVars[1], "% Var Expl"),
    ylim = range(pca$x[, 1]),
    outline = FALSE
)
points(pca$x[, 1] ~ jitter(as.numeric(rse_gene$COI), amount = 0.15),
    pch = 21,
    bg = rse_gene$COI
)
dev.off()

pdf(file.path(dir_plots, "PC1_vs_totalAssignedGene.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
plot(
    pca$x[, 1] ~ colData(rse_gene)[, nn[1]],
    xlab = nn[1],
    ylab = paste0("PC1: ", pcaVars[1], "% Var Expl"),
    pch = 21,
    bg = rse_gene$COI
)
for (i in seq(along = gIndexes)) {
    ii <- gIndexes[[i]]
    abline(lm(pca$x[, 1] ~ colData(rse_gene)[, nn[1]], subset = ii),
        lwd = 4,
        col = i
    )
}

dev.off()

pdf(file.path(dir_plots, "totalAssignedGene_vs_COI.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
boxplot(
    colData(rse_gene)[, nn[1]] ~ rse_gene$COI,
    ylab = nn[1],
    xlab = "",
    outline = FALSE,
    ylim = range(colData(rse_gene)[, nn[1]])
)
points(
    colData(rse_gene)[, nn[1]] ~
        jitter(as.numeric(rse_gene$COI), amount = 0.15),
    pch = 21,
    bg = rse_gene$COI
)
dev.off()

## PC2
pdf(file.path(dir_plots, "PC2_vs_totalAssignedGene.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
plot(
    pca$x[, 2] ~ colData(rse_gene)[, nn[1]],
    xlab = nn[1],
    ylab = paste0("PC2: ", pcaVars[2], "% Var Expl"),
    pch = 21,
    bg = rse_gene$COI
)
for (i in seq(along = gIndexes)) {
    ii <- gIndexes[[i]]
    abline(lm(pca$x[, 2] ~ colData(rse_gene)[, nn[1]], subset = ii),
        lwd = 4,
        col = i
    )
}

dev.off()


pdf(file.path(dir_plots, "PC2_vs_COI.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
boxplot(
    pca$x[, 2] ~ rse_gene$COI,
    xlab = "",
    ylab = paste0("PC2: ", pcaVars[2], "% Var Expl"),
    ylim = range(pca$x[, 2]),
    outline = FALSE
)
points(pca$x[, 2] ~ jitter(as.numeric(rse_gene$COI), amount = 0.15),
    pch = 21,
    bg = rse_gene$COI
)
dev.off()

### PC3
pdf(file.path(dir_plots, "PC3_vs_COI.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
boxplot(
    pca$x[, 3] ~ rse_gene$COI,
    xlab = "",
    ylab = paste0("PC3: ", pcaVars[3], "% Var Expl"),
    ylim = range(pca$x[, 3]),
    outline = FALSE
)
points(pca$x[, 3] ~ jitter(as.numeric(rse_gene$COI), amount = 0.15),
    pch = 21,
    bg = rse_gene$COI
)
dev.off()

pdf(file.path(dir_plots, "PC3_vs_totalAssignedGene.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
plot(
    pca$x[, 3] ~ colData(rse_gene)[, nn[1]],
    xlab = nn[1],
    ylab = paste0("PC3: ", pcaVars[3], "% Var Expl"),
    pch = 21,
    bg = rse_gene$COI
)
for (i in seq(along = gIndexes)) {
    ii <- gIndexes[[i]]
    abline(lm(pca$x[, 3] ~ colData(rse_gene)[, nn[1]], subset = ii),
        lwd = 4,
        col = i
    )
}
dev.off()


pdf(file.path(dir_plots, "mitoRate_vs_COI.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 1.7,
    cex.lab = 1.7
)
palette(brewer.pal(4, "Dark2"))
boxplot(
    rse_gene$mitoRate ~ rse_gene$COI,
    xlab = "",
    ylab = paste0("mitoRate"),
    ylim = range(rse_gene$mitoRate),
    outline = FALSE
)
points(rse_gene$mitoRate ~ jitter(as.numeric(rse_gene$COI), amount = 0.15),
    pch = 21,
    bg = rse_gene$COI
)
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
#  Biobase              * 2.60.0    2023-04-25 [1] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [1] Bioconductor
#  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  brio                   1.1.3     2021-11-30 [1] CRAN (R 4.3.0)
#  cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  callr                  3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  data.table             1.14.8    2023-02-17 [1] CRAN (R 4.3.0)
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
#  IRanges              * 2.34.1    2023-07-02 [1] Bioconductor
#  jaffelab             * 0.99.32   2023-05-06 [1] Github (LieberInstitute/jaffelab@7b7afe3)
#  later                  1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
#  lattice              * 0.21-8    2023-04-05 [1] CRAN (R 4.3.1)
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
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
#  nlme                   3.1-163   2023-08-09 [1] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild               1.4.2     2023-06-26 [1] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload                1.3.2.1   2023-07-08 [1] CRAN (R 4.3.0)
#  prettyunits            1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
#  processx               3.8.2     2023-06-30 [1] CRAN (R 4.3.0)
#  profvis                0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
#  promises               1.2.1     2023-08-10 [1] CRAN (R 4.3.0)
#  prompt                 1.0.2     2023-08-31 [1] CRAN (R 4.3.0)
#  ps                     1.7.5     2023-04-18 [1] CRAN (R 4.3.0)
#  purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
#  RColorBrewer         * 1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
#  RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
#  readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
#  remotes                2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
#  rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
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
#  testthat             * 3.1.10    2023-07-06 [1] CRAN (R 4.3.0)
#  tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
#  timechange             0.2.0     2023-01-11 [1] CRAN (R 4.3.0)
#  urlchecker             1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
#  usethis              * 2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
#  vctrs                  0.6.3     2023-06-14 [1] CRAN (R 4.3.0)
#  xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-04-25 [1] Bioconductor
#  zlibbioc               1.46.0    2023-04-25 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
