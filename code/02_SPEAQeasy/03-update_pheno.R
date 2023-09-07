library("SummarizedExperiment")
library("here")
library("readr")
library("dplyr")
library("sessioninfo")

## Read in the individual data
individual <- read_csv(here("raw-data/psychENCODE/v2023-07-14/metadata/LBP_individual_metadata.csv"))

## Some columns are completely empty
colSums(is.na(individual)) == nrow(individual)
# individualID individualIdSource            species                sex               race          ethnicity     yearsEducation
#        FALSE              FALSE              FALSE              FALSE              FALSE              FALSE               TRUE
#     ageDeath         causeDeath        mannerDeath       apoeGenotype                pmi                 pH        brainWeight
#        FALSE               TRUE               TRUE               TRUE              FALSE               TRUE               TRUE
#    diagnosis  diagnosisCriteria              CERAD              Braak
#        FALSE               TRUE               TRUE               TRUE
individual <- individual[, colSums(is.na(individual)) != nrow(individual)]
individual
# # A tibble: 411 × 9
#    individualID       individualIdSource species sex    race  ethnicity ageDeath pmi    diagnosis
#    <chr>              <chr>              <chr>   <chr>  <chr> <chr>     <chr>    <time> <chr>
#  1 HBDK_17_002_FC     BEB_Miller         Human   male   White NHL       79          NA  Parkinson's disease
#  2 HBGM_17_001_FC     BEB_Miller         Human   male   White NHL       89+         NA  Parkinson's disease
#  3 HBGQ_17_001_FC     BEB_Miller         Human   female White NHL       89+         NA  Parkinson's disease
#  4 HBHE_17_001_FC     BEB_Miller         Human   female White NHL       81          NA  Parkinson's disease
#  5 HBHJ_17_001_FC     BEB_Miller         Human   female White NHL       79          NA  Parkinson's disease
#  6 HBHV_17_001_FC     BEB_Miller         Human   female White NHL       79          NA  Parkinson's disease
#  7 HBIE_17_001_FC     BEB_Miller         Human   male   White NHL       80          NA  Parkinson's disease
#  8 HBJE_17_001_FC     BEB_Miller         Human   male   White NHL       80          NA  Parkinson's disease
#  9 Hct15HAG_17_001_FC BEB_Miller         Human   female White NHL       86          NA  control
# 10 Hct15HAL_17_003_FC BEB_Miller         Human   male   White NHL       61          NA  control
# # ℹ 401 more rows
# # ℹ Use `print(n = ...)` to see more rows

## Let's update some of them
individual <- mutate(individual,
    individualIdSource = factor(individualIdSource),
    species = factor(species),
    sex = factor(sex),
    race = factor(race),
    ethnicity = factor(ethnicity),
    ageDeath_num = as.numeric(gsub("\\+", "", ageDeath)),
    diagnosis = factor(diagnosis)
)
individual <- mutate(individual,
    diagnosis = relevel(diagnosis, "control")
)
individual
# # A tibble: 411 × 10
#    individualID       individualIdSource species sex    race  ethnicity ageDeath pmi    diagnosis           ageDeath_num
#    <chr>              <fct>              <fct>   <fct>  <fct> <fct>     <chr>    <time> <fct>                      <dbl>
#  1 HBDK_17_002_FC     BEB_Miller         Human   male   White NHL       79          NA  Parkinson's disease           79
#  2 HBGM_17_001_FC     BEB_Miller         Human   male   White NHL       89+         NA  Parkinson's disease           89
#  3 HBGQ_17_001_FC     BEB_Miller         Human   female White NHL       89+         NA  Parkinson's disease           89
#  4 HBHE_17_001_FC     BEB_Miller         Human   female White NHL       81          NA  Parkinson's disease           81
#  5 HBHJ_17_001_FC     BEB_Miller         Human   female White NHL       79          NA  Parkinson's disease           79
#  6 HBHV_17_001_FC     BEB_Miller         Human   female White NHL       79          NA  Parkinson's disease           79
#  7 HBIE_17_001_FC     BEB_Miller         Human   male   White NHL       80          NA  Parkinson's disease           80
#  8 HBJE_17_001_FC     BEB_Miller         Human   male   White NHL       80          NA  Parkinson's disease           80
#  9 Hct15HAG_17_001_FC BEB_Miller         Human   female White NHL       86          NA  control                       86
# 10 Hct15HAL_17_003_FC BEB_Miller         Human   male   White NHL       61          NA  control                       61
# # ℹ 401 more rows
# # ℹ Use `print(n = ...)` to see more rows
summary(individual)
# individualID          individualIdSource  species        sex                                    race
# Length:411         BEB_Miller  : 13      Human:411   female:156   American Indian or Alaska Native:  1
# Class :character   HBTRC       :104                  male  :255   Asian                           :  7
# Mode  :character   MSSM_Charney:168                               Black or African American       :  4
#                    NYBB        :126                               White                           :199
#                                                                   NA's                            :200
#
#
#     ethnicity     ageDeath             pmi                                   diagnosis    ageDeath_num
# NHL      :158   Length:411         Length:411        control                      : 96   Min.   :24.0
# HL       :  9   Class :character   Class1:hms        Alzheimer Disease            : 14   1st Qu.:62.0
# English  :  2   Mode  :character   Class2:difftime   dystonia                     : 13   Median :70.0
# Caucasion:  1                      Mode  :numeric    essential tremor             : 17   Mean   :68.7
# German   :  1                                        Major Depressive Disorder    :  2   3rd Qu.:77.0
# (Other)  :  4                                        obsessive compulsive disorder:  5   Max.   :89.0
# NA's     :236                                        Parkinson's disease          :264   NA's   :2

## Read in the biospecimen info which we can use to match the specimenID to the individualID
biospecimen <- read_csv(here("raw-data/psychENCODE/v2023-07-14/metadata/LBP_biospecimen_metadata.csv"))

## Some columns are completely empty
colSums(is.na(biospecimen)) == nrow(biospecimen)
# individualID        specimenID  specimenIdSource             organ            tissue      BrodmannArea
#        FALSE             FALSE             FALSE             FALSE             FALSE              TRUE
# sampleStatus      tissueWeight      tissueVolume nucleicAcidSource          cellType      fastingState
#        FALSE              TRUE              TRUE             FALSE              TRUE              TRUE
# isPostMortem       samplingAge  samplingAgeUnits       visitNumber             assay
#        FALSE             FALSE             FALSE              TRUE             FALSE

biospecimen <- biospecimen[, colSums(is.na(biospecimen)) != nrow(biospecimen)]
biospecimen <- mutate(biospecimen,
    specimenIdSource = factor(specimenIdSource),
    organ = factor(organ),
    tissue = factor(tissue),
    sampleStatus = factor(sampleStatus),
    nucleicAcidSource = factor(nucleicAcidSource),
    ageDeath_num = as.numeric(gsub("\\+", "", samplingAge)),
    samplingAgeUnits = factor(samplingAgeUnits),
    assay = factor(assay),
    samplingAge = NULL
)
biospecimen
# # A tibble: 516 × 11
#    individualID       specimenID       specimenIdSource organ tissue    sampleStatus nucleicAcidSource isPostMortem samplingAgeUnits assay ageDeath_num
#    <chr>              <chr>            <fct>            <fct> <fct>     <fct>        <fct>             <lgl>        <fct>            <fct>        <dbl>
#  1 HBDK_17_002_FC     LBPSEMA4BRAIN518 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           79
#  2 HBGM_17_001_FC     LBPSEMA4BRAIN036 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           89
#  3 HBGQ_17_001_FC     LBPSEMA4BRAIN273 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           89
#  4 HBHE_17_001_FC     LBPSEMA4BRAIN585 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           81
#  5 HBHJ_17_001_FC     LBPSEMA4BRAIN535 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           79
#  6 HBHV_17_001_FC     LBPSEMA4BRAIN392 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           79
#  7 HBIE_17_001_FC     LBPSEMA4BRAIN522 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           80
#  8 HBJE_17_001_FC     LBPSEMA4BRAIN358 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           80
#  9 Hct15HAG_17_001_FC LBPSEMA4BRAIN723 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           86
# 10 Hct15HAL_17_003_FC LBPSEMA4BRAIN220 BEB_Miller       brain dorsolat… frozen       bulk cell         TRUE         years            rnaS…           61
# # ℹ 506 more rows
# # ℹ Use `print(n = ...)` to see more rows

summary(biospecimen)
# individualID        specimenID            specimenIdSource   organ                                tissue    sampleStatus
# Length:516         Length:516         BEB_Miller  : 13     brain:516   dorsolateral prefrontal cortex: 13   frozen:516
# Class :character   Class :character   HBTRC       :104                 left cerebral hemisphere      :269
# Mode  :character   Mode  :character   MSSM_Charney:273                 right cerebral hemisphere     :234
#                                       NYBB        :126
#
#
#
# nucleicAcidSource isPostMortem    samplingAgeUnits    assay      ageDeath_num
# bulk cell:516     Mode :logical   years:516        rnaSeq:516   Min.   :24.00
#                   FALSE:273                                     1st Qu.:60.25
#                   TRUE :243                                     Median :68.00
#                                                                 Mean   :67.10
#                                                                 3rd Qu.:76.00
#                                                                 Max.   :89.00
#                                                                 NA's   :2

synapse <- left_join(biospecimen, individual, by = "individualID", suffix = c(".biospecimen", ".individual"))

## There's one mistmatch by age. The individual table says it's 75 years old
## whereas the biospecimen table says it was 74 years old.
synapse[ which(synapse$ageDeath_num.biospecimen - synapse$ageDeath_num.individual == -1), c("individualID", "specimenID", "ageDeath_num.biospecimen", "ageDeath_num.individual")]
# # A tibble: 1 × 4
#   individualID specimenID       ageDeath_num.biospecimen ageDeath_num.individual
#   <chr>        <chr>                               <dbl>                   <dbl>
# 1 PT-0058      LBPSEMA4BRAIN553                       74                      75

## Read in the RNAseq table
rnaseq <- read_csv(here("raw-data/psychENCODE/v2023-07-14/metadata/LBP_assay_RNAseq_metadata.csv"))

## Keep only the columns with some data
rnaseq <- rnaseq[, colSums(is.na(rnaseq)) != nrow(rnaseq)]

## Update some columns
rnaseq <- mutate(rnaseq,
    platform = factor(platform),
    rnaBatch = factor(rnaBatch),
    libraryBatch = factor(libraryBatch),
    sequencingBatch = factor(sequencingBatch),
    libraryPrep = factor(libraryPrep),
    libraryPreparationMethod = factor(libraryPreparationMethod),
    readStrandOrigin = factor(readStrandOrigin),
    runType = factor(runType),
    assay = factor(assay)
)
rnaseq
# # A tibble: 516 × 13
#    specimenID       platform              RIN rnaBatch         libraryBatch sequencingBatch libraryPrep   libraryPreparationMethod isStranded readStrandOrigin runType   readLength assay
#    <chr>            <fct>               <dbl> <fct>            <fct>        <fct>           <fct>         <fct>                    <lgl>      <fct>            <fct>          <dbl> <fct>
#  1 LBPSEMA4BRAIN518 IlluminaNovaseq6000   7.1 LBPSEMA4_PLATE_1 db3          sb5             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  2 LBPSEMA4BRAIN036 IlluminaNovaseq6000   7.5 LBPSEMA4_PLATE_5 db20         sb2             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  3 LBPSEMA4BRAIN273 IlluminaNovaseq6000   6.4 LBPSEMA4_PLATE_5 db9          sb7             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  4 LBPSEMA4BRAIN585 IlluminaNovaseq6000   7.1 LBPSEMA4_PLATE_4 db11         sb3             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  5 LBPSEMA4BRAIN535 IlluminaNovaseq6000   6.2 LBPSEMA4_PLATE_7 db22         sb6             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  6 LBPSEMA4BRAIN392 IlluminaNovaseq6000   5.8 LBPSEMA4_PLATE_7 db14         sb4             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  7 LBPSEMA4BRAIN522 IlluminaNovaseq6000   6.9 LBPSEMA4_PLATE_6 db18         sb8             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  8 LBPSEMA4BRAIN358 IlluminaNovaseq6000   6.9 LBPSEMA4_PLATE_4 db14         sb4             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
#  9 LBPSEMA4BRAIN723 IlluminaNovaseq6000   7   LBPSEMA4_PLATE_6 db8          sb7             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
# 10 LBPSEMA4BRAIN220 IlluminaNovaseq6000   8.7 LBPSEMA4_PLATE_7 db22         sb6             rRNAdepletion TruSeq                   TRUE       reverse          pairedEnd        100 rnaSeq
# # ℹ 506 more rows
# # ℹ Use `print(n = ...)` to see more rows
summary(rnaseq)
# specimenID                       platform        RIN                    rnaBatch    libraryBatch sequencingBatch        libraryPrep
# Length:516         IlluminaNovaseq6000:516   Min.   :4.100   LBPSEMA4_PLATE_4: 64   db17   : 23   sb5    : 65     rRNAdepletion:516
# Class :character                             1st Qu.:6.700   LBPSEMA4_PLATE_6: 64   db10   : 22   sb6    : 65
# Mode  :character                             Median :7.300   LBPSEMA4_PLATE_7: 63   db13   : 22   sb2    : 64
#                                              Mean   :7.296   LBPSEMA4_PLATE_1: 61   db2    : 22   sb4    : 63
#                                              3rd Qu.:7.825   LBPSEMA4_PLATE_2: 61   db21   : 22   sb1    : 62
#                                              Max.   :9.700   LBPSEMA4_PLATE_5: 61   db23   : 22   sb3    : 62
#                                                              (Other)         :142   (Other):383   (Other):135
# libraryPreparationMethod isStranded     readStrandOrigin      runType      readLength     assay
# TruSeq:516               Mode:logical   reverse:516      pairedEnd:516   Min.   :100   rnaSeq:516
#                          TRUE:516                                        1st Qu.:100
#                                                                          Median :100
#                                                                          Mean   :100
#                                                                          3rd Qu.:100
#                                                                          Max.   :100

synapse <- left_join(synapse, rnaseq, by = "specimenID")
stopifnot(with(synapse, identical(assay.x, assay.y)))
synapse <- mutate(synapse,
    assay = assay.x,
    assay.x = NULL,
    assay.y = NULL
)
summary(synapse)
# individualID        specimenID            specimenIdSource   organ                                tissue    sampleStatus
# Length:516         Length:516         BEB_Miller  : 13     brain:516   dorsolateral prefrontal cortex: 13   frozen:516
# Class :character   Class :character   HBTRC       :104                 left cerebral hemisphere      :269
# Mode  :character   Mode  :character   MSSM_Charney:273                 right cerebral hemisphere     :234
#                                       NYBB        :126
#
#
#
# nucleicAcidSource isPostMortem    samplingAgeUnits ageDeath_num.biospecimen    individualIdSource  species        sex
# bulk cell:516     Mode :logical   years:516        Min.   :24.00            BEB_Miller  : 13      Human:516   female:197
#                   FALSE:273                        1st Qu.:60.25            HBTRC       :104                  male  :319
#                   TRUE :243                        Median :68.00            MSSM_Charney:273
#                                                    Mean   :67.10            NYBB        :126
#                                                    3rd Qu.:76.00
#                                                    Max.   :89.00
#                                                    NA's   :2
#                               race         ethnicity     ageDeath             pmi                                   diagnosis
# American Indian or Alaska Native:  1   NHL      :250   Length:516         Length:516        control                      : 96
# Asian                           : 13   HL       : 16   Class :character   Class1:hms        Alzheimer Disease            : 14
# Black or African American       :  4   English  :  2   Mode  :character   Class2:difftime   dystonia                     : 20
# White                           :292   Indian   :  2                      Mode  :numeric    essential tremor             : 21
# NA's                            :206   Caucasion:  1                                        Major Depressive Disorder    :  4
#                                        (Other)  :  4                                        obsessive compulsive disorder:  9
#                                        NA's     :241                                        Parkinson's disease          :352
# ageDeath_num.individual                platform        RIN                    rnaBatch    libraryBatch sequencingBatch
# Min.   :24.00           IlluminaNovaseq6000:516   Min.   :4.100   LBPSEMA4_PLATE_4: 64   db17   : 23   sb5    : 65
# 1st Qu.:60.25                                     1st Qu.:6.700   LBPSEMA4_PLATE_6: 64   db10   : 22   sb6    : 65
# Median :68.00                                     Median :7.300   LBPSEMA4_PLATE_7: 63   db13   : 22   sb2    : 64
# Mean   :67.10                                     Mean   :7.296   LBPSEMA4_PLATE_1: 61   db2    : 22   sb4    : 63
# 3rd Qu.:76.00                                     3rd Qu.:7.825   LBPSEMA4_PLATE_2: 61   db21   : 22   sb1    : 62
# Max.   :89.00                                     Max.   :9.700   LBPSEMA4_PLATE_5: 61   db23   : 22   sb3    : 62
# NA's   :2                                                         (Other)         :142   (Other):383   (Other):135
#        libraryPrep  libraryPreparationMethod isStranded     readStrandOrigin      runType      readLength     assay
# rRNAdepletion:516   TruSeq:516               Mode:logical   reverse:516      pairedEnd:516   Min.   :100   rnaSeq:516
#                                              TRUE:516                                        1st Qu.:100
#                                                                                              Median :100
#                                                                                              Mean   :100
#                                                                                              3rd Qu.:100
#                                                                                              Max.   :100

## Load the rse_gene object
load(here("processed-data", "02_SPEAQeasy", "pipeline_output", "count_objects", "rse_gene_living_brain_reanalysis_n516.Rdata"), verbose = TRUE)

## Extract the phenotype information
pheno <- colData(rse_gene)

## Merge the SPEAQeasy outputs with the information available from Synapse
pheno_with_synapse <- DataFrame(left_join(as.data.frame(pheno), synapse, by = join_by("SAMPLE_ID" == "specimenID")))
rownames(pheno_with_synapse) <- rownames(pheno)

## Replace the colData() info
colData(rse_gene) <- pheno_with_synapse

## Save for later
write_csv(as.data.frame(pheno_with_synapse), here("processed-data/02_SPEAQeasy/", "pheno_with_synapse.csv"))
saveRDS(pheno_with_synapse, file = here("processed-data/02_SPEAQeasy/", "pheno_with_synapse.Rds"))
saveRDS(rse_gene, file = here("processed-data/02_SPEAQeasy/", "rse_gene_living_brain_reanalysis_n516.Rds"))

## Repeat for the other feature types
rm(rse_gene)
load(here("processed-data", "02_SPEAQeasy", "pipeline_output", "count_objects", "rse_exon_living_brain_reanalysis_n516.Rdata"), verbose = TRUE)
colData(rse_exon) <- pheno_with_synapse
saveRDS(rse_exon, file = here("processed-data/02_SPEAQeasy/", "rse_exon_living_brain_reanalysis_n516.Rds"))

rm(rse_exon)
load(here("processed-data", "02_SPEAQeasy", "pipeline_output", "count_objects", "rse_jx_living_brain_reanalysis_n516.Rdata"), verbose = TRUE)
colData(rse_jx) <- pheno_with_synapse
saveRDS(rse_jx, file = here("processed-data/02_SPEAQeasy/", "rse_jx_living_brain_reanalysis_n516.Rds"))

rm(rse_jx)
load(here("processed-data", "02_SPEAQeasy", "pipeline_output", "count_objects", "rse_tx_living_brain_reanalysis_n516.Rdata"), verbose = TRUE)
colData(rse_tx) <- pheno_with_synapse
saveRDS(rse_tx, file = here("processed-data/02_SPEAQeasy/", "rse_tx_living_brain_reanalysis_n516.Rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# sgejobs::job_single("03-update_pheno", create_shell = TRUE, queue = "caracol", memory = "100G", command = "Rscript 03-update_pheno.R")
