library("here")
library("recount3")
library("sessioninfo")

## Create output directories
dir_rdata <-
    here("processed-data", "04_check_degradation", "recount3_salmon")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)


## Locate Salmon output from recount3 for SRP108559
salmon_files <-
    list.files(
        "/dcl02/lieber/ajaffe/trash_jhpce_org/recount-pump/human_tranche_backups/sra_human_v3_9/59/SRP108559",
        pattern = "sra.salmon.tsv.zst$",
        full.names = TRUE,
        recursive = TRUE
    )
stopifnot(length(salmon_files) == 40)

## Copy to this project
file.copy(salmon_files, file.path(dir_rdata, basename(salmon_files)))

## Uncompress the local files
local_salmon_zst <-
    list.files(dir_rdata, pattern = "sra.salmon.tsv.zst$", full.names = TRUE)
sapply(paste("unzstd", local_salmon_zst), system)

## Read the local files
local_salmon <-
    list.files(dir_rdata, pattern = "sra.salmon.tsv$", full.names = TRUE)
salmon_list <- lapply(local_salmon, read.table, header = TRUE)
names(salmon_list) <- gsub("\\!.*", "", basename(local_salmon))

## Delete uncompressed files
unlink(local_salmon)

## Download the metadata from recount3 for SRP108559
hp <- available_projects()
rse_degrade <- create_rse(hp[hp$project == "SRP108559", ])
rse_degrade <- expand_sra_attributes(rse_degrade)

## Re-order the salmon output to match the recount3 order
stopifnot(all(names(salmon_list) %in% colnames(rse_degrade)))
salmon_list <- salmon_list[colnames(rse_degrade)]

## Extract Assays
TPM <- do.call(cbind, lapply(salmon_list, "[[", "TPM"))
EffectiveLength <-
    do.call(cbind, lapply(salmon_list, "[[", "EffectiveLength"))
NumReads <- do.call(cbind, lapply(salmon_list, "[[", "NumReads"))

## Check that the are all the same length (aka, same order)
Length <- do.call(cbind, lapply(salmon_list, "[[", "Length"))
stopifnot(all(rowMeans(Length) == Length[, 1]))

## Build an RSE for later use
rse_degrade_tx <- SummarizedExperiment(
    assays = SimpleList(
        TPM = TPM,
        EffectiveLength = EffectiveLength,
        NumReads = NumReads
    ),
    rowData = DataFrame(
        gencode_id = salmon_list[[1]]$Name,
        ensembl_id = gsub("\\..*", "", salmon_list[[1]]$Name),
        Length = Length[, 1]
    ),
    colData = colData(rse_degrade),
    metadata = metadata(rse_degrade)
)
rownames(rse_degrade_tx) <- rowData(rse_degrade_tx)$gencode_id

## Save for later
saveRDS(
    rse_degrade_tx,
    file = here(
        "processed-data",
        "04_check_degradation",
        "rse_degrade_tx.Rds"
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 Patched (2023-09-08 r85106)
#  os       CentOS Linux 7 (Core)
#  system   x86                 1.10.0    2023-04-25 [2] Bioconductor
#  BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
#  Biostrings             2.68.1    2023-05-16 [2] Bioconductor
#  bit                    4.0.5     2022-11-15 [2] CRAN (R 4.3.0)
#  bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.3.0)
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.0)
#  blob                   1.2.4     2023-03-17 [2] CRAN (R 4.3.0)
#  cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.0)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.0)
#  curl                   5.0.2     2023-08-14 [2] CRAN (R 4.3.1)
#  data.table             1.14.8    2023-02-17 [2] CRAN (R 4.3.0)
#  DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.3.0)
#  dbplyr                 2.3.3     2023-07-07 [2] CRAN (R 4.3.1)
#  DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
#  digest                 0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
#  dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
#  fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.0)
#  filelock               1.0.2     2018-10-05 [2] CRAN (R 4.3.0)
#  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.2    2023-08-25 [2] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-04-11 [2] Bioconductor
#  GenomicAlignments      1.36.0    2023-04-25 [2] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
#  ggplot2                3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
#  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.0)
#  gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.0)
#  htmltools              0.5.6     2023-08-10 [2] CRAN (R 4.3.1)
#  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.3.0)
#  httpuv                 1.6.11    2023-05-11 [2] CRAN (R 4.3.0)
#  httr                   1.4.7     2023-08-15 [2] CRAN (R 4.3.1)
#  IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
#  jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
#  later                  1.3.1     2023-05-02 [2] CRAN (R 4.3.0)
#  lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
#  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.0)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.0)
#  Matrix                 1.6-1     2023-08-14 [3] CRAN (R 4.3.1)
#  MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.0)
#  memoise                2.0.1     2021-11-26 [2] CRAN (R 4.3.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.0)
#  png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.0)
#  promises               1.2.1     2023-08-10 [2] CRAN (R 4.3.1)
#  purrr                  1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
#  R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.0)
#  R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.0)
#  R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
#  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.0)
#  recount3             * 1.10.2    2023-05-07 [2] Bioconductor
#  restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.0)
#  rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.0)
#  rmote                  0.3.4     2023-05-06 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.0)
#  Rsamtools              2.16.0    2023-04-25 [2] Bioconductor
#  RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.3.0)
#  rtracklayer            1.60.1    2023-08-15 [2] Bioconductor
#  S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
#  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.0)
#  servr                  0.27      2023-05-02 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.0)
#  SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.0)
#  vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
#  withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.0)
#  xfun                   0.40      2023-08-09 [2] CRAN (R 4.3.1)
#  XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-04-25 [2] Bioconductor
#  yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.0)
#  zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
#
#  [1] /users/lcollado/R/4.3
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
