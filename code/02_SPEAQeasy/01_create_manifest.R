library("here")
library("sessioninfo")

dir_rdata <- here("processed-data", "02_SPEAQeasy")
dir.create(dir_rdata, recursive = TRUE, showWarnings = FALSE)

## Read in information from the Synapse download query
synapse <-
    read.csv(
        here(
            "raw-data/psychENCODE/v2023-07-14/RNAseq/SYNAPSE_TABLE_QUERY_127308191.csv"
        )
    )
## Subset to just RNAseq data
synapse <- subset(synapse, fileFormat == "fastq")
dim(synapse)
# [1] 7838   46

## Check that individualID and specimenID do not overlap
table(sapply(with(synapse, split(
    individualID, specimenID
)), function(x) {
    length(unique(x))
}))
#   1
# 516

## Check that all fastq.gz files were downloaded
synapse$JHPCE_path <-
    file.path(
        here("raw-data/psychENCODE/v2023-07-14/RNAseq"),
        synapse$name
    )
synapse$JHPCE_file_exists <- file.exists(synapse$JHPCE_path)
table(synapse$JHPCE_file_exists)
# FALSE  TRUE
#   227  7611
synapse_missing <- subset(synapse, !JHPCE_file_exists)
synapse_present <- subset(synapse, JHPCE_file_exists)

unique(synapse_missing$specimenID)
#  [1] "LBPSEMA4BRAIN639" "LBPSEMA4BRAIN730" "LBPSEMA4BRAIN176" "LBPSEMA4BRAIN285"
#  [5] "LBPSEMA4BRAIN469" "LBPSEMA4BRAIN043" "LBPSEMA4BRAIN141" "LBPSEMA4BRAIN085"
#  [9] "LBPSEMA4BRAIN607" "LBPSEMA4BRAIN401" "LBPSEMA4BRAIN348" "LBPSEMA4BRAIN478"
# [13] "LBPSEMA4BRAIN402" "LBPSEMA4BRAIN778" "LBPSEMA4BRAIN623"
unique(synapse_missing$specimenID) %in% synapse_present$specimenID
#  [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [13] FALSE FALSE FALSE

## Official syntax from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520
# synapse get -q "SELECT * FROM syn11346063.39 WHERE ( \"study\" HAS ( 'LBP' ) ) AND ( \"specimenID\" LIKE '%LBPSEMA4BRAIN639%' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN730%' )"
paste(paste0(
    "\"specimenID\" LIKE '%",
    unique(synapse_missing$specimenID),
    "' OR"
), collapse = " ")
# "\"specimenID\" LIKE '%LBPSEMA4BRAIN639' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN730' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN176' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN285' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN469' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN043' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN141' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN085' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN607' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN401' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN348' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN478' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN402' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN778' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN623' OR"


## Read in the data from the missing samples
synapse_v2 <-
    read.csv(
        here(
            "raw-data/psychENCODE/v2023-07-31/RNAseq/SYNAPSE_TABLE_QUERY_127507067.csv"
        )
    )
dim(synapse_v2)
# [1] 228  46
synapse_v2$JHPCE_path <-
    file.path(
        here("raw-data/psychENCODE/v2023-07-31/RNAseq"),
        synapse_v2$name
    )
synapse_v2$JHPCE_file_exists <- file.exists(synapse_v2$JHPCE_path)
table(synapse_v2$JHPCE_file_exists)
# TRUE
#  228

## Check that the new data was indeed the one missing before
table(subset(synapse, name %in% synapse_v2$name)$JHPCE_file_exists)
# FALSE  TRUE
#   227     1

## Merge all the data together
synapse_full <- rbind(
    subset(synapse, !name %in% synapse_v2$name),
    synapse_v2
)
stopifnot(all(synapse_full$JHPCE_file_exists))

## Export for later use
write.csv(synapse_full, file = file.path(dir_rdata, "synapse_full.csv"))

## Get read prefix
synapse_full$name_prefix <-
    gsub("_[[:digit:]]+.fastq.gz", "", synapse_full$name)
synapse_full$name_run <-
    gsub("_R[[:digit:]]$", "", synapse_full$name_prefix)
synapse_full$name_read_pair <-
    gsub(".*_", "", synapse_full$name_prefix)
table(synapse_full$name_read_pair)
#   R1   R2
# 3919 3919

## Split by read pair
synapse_full_r1 <- subset(synapse_full, name_read_pair == "R1")
synapse_full_r2 <- subset(synapse_full, name_read_pair == "R2")
## Check that they are in the same order
m <- match(synapse_full_r1$name_run, synapse_full_r2$name_run)
stopifnot(identical(m, seq_len(length(m))))

## Build the manifest file for SPEAQeasy
manifest <- data.frame(
    read1 = synapse_full_r1$JHPCE_path,
    read1_md5sum = rep(0, nrow(synapse_full_r1)),
    read2 = synapse_full_r2$JHPCE_path,
    read2_md5sum = rep(0, nrow(synapse_full_r1)),
    sample_name = synapse_full_r1$specimenID
)
dim(manifest)
# [1] 3919    5
length(unique(manifest$sample_name))
# [1] 516

## Export the manifest file for SPEAQeasy
write.table(
    manifest,
    file = file.path(dir_rdata, "samples.manifest"),
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 Patched (2023-08-09 r84918)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-08-10
#  pandoc   3.1.1 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/bin/pandoc
#
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date (UTC) lib source
#  cli           3.6.1   2023-03-23 [2] CRAN (R 4.3.0)
#  colorout      1.2-2   2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace    2.1-0   2023-01-23 [2] CRAN (R 4.3.0)
#  digest        0.6.33  2023-07-07 [2] CRAN (R 4.3.1)
#  dplyr         1.1.2   2023-04-20 [2] CRAN (R 4.3.0)
#  fansi         1.0.4   2023-01-22 [2] CRAN (R 4.3.0)
#  fastmap       1.1.1   2023-02-24 [2] CRAN (R 4.3.0)
#  generics      0.1.3   2022-07-05 [2] CRAN (R 4.3.0)
#  ggplot2       3.4.2   2023-04-03 [2] CRAN (R 4.3.0)
#  glue          1.6.2   2022-02-24 [2] CRAN (R 4.3.0)
#  gtable        0.3.3   2023-03-21 [2] CRAN (R 4.3.0)
#  here        * 1.0.1   2020-12-13 [2] CRAN (R 4.3.0)
#  htmltools     0.5.5   2023-03-23 [2] CRAN (R 4.3.0)
#  htmlwidgets   1.6.2   2023-03-17 [2] CRAN (R 4.3.0)
#  httpuv        1.6.11  2023-05-11 [2] CRAN (R 4.3.0)
#  jsonlite      1.8.7   2023-06-29 [2] CRAN (R 4.3.1)
#  later         1.3.1   2023-05-02 [2] CRAN (R 4.3.0)
#  lattice       0.21-8  2023-04-05 [3] CRAN (R 4.3.1)
#  lifecycle     1.0.3   2022-10-07 [2] CRAN (R 4.3.0)
#  magrittr      2.0.3   2022-03-30 [2] CRAN (R 4.3.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 4.3.0)
#  pillar        1.9.0   2023-03-22 [2] CRAN (R 4.3.0)
#  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.3.0)
#  png           0.1-8   2022-11-29 [2] CRAN (R 4.3.0)
#  promises      1.2.0.1 2021-02-11 [2] CRAN (R 4.3.0)
#  R6            2.5.1   2021-08-19 [2] CRAN (R 4.3.0)
#  Rcpp          1.0.11  2023-07-06 [2] CRAN (R 4.3.1)
#  rlang         1.1.1   2023-04-28 [2] CRAN (R 4.3.0)
#  rmote         0.3.4   2023-05-06 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot     2.0.3   2022-04-02 [2] CRAN (R 4.3.0)
#  scales        1.2.1   2022-08-20 [2] CRAN (R 4.3.0)
#  servr         0.27    2023-05-02 [1] CRAN (R 4.3.0)
#  sessioninfo * 1.2.2   2021-12-06 [2] CRAN (R 4.3.0)
#  tibble        3.2.1   2023-03-20 [2] CRAN (R 4.3.0)
#  tidyselect    1.2.0   2022-10-10 [2] CRAN (R 4.3.0)
#  utf8          1.2.3   2023-01-31 [2] CRAN (R 4.3.0)
#  vctrs         0.6.3   2023-06-14 [2] CRAN (R 4.3.1)
#  xfun          0.39    2023-04-20 [2] CRAN (R 4.3.0)
#
#  [1] /users/lcollado/R/4.3
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/library
#
# ────────────────────────────────────────────────────────────────────────────────────────────────────────
