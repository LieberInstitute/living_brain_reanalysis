library("qsvaR")
library("SummarizedExperiment")
library("parallel")
library("jaffelab")
library("limma")
library("rtracklayer")
library("RColorBrewer")
library("here")
library("sessioninfo")

## For reproducing the jitter output
set.seed(20230907)

## Create output directories
dir_plots <- here("plots", "04_check_degradation")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir_rdata <- here("processed-data", "04_check_degradation")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load SPEAQeasy gene level data
rse_tx_full <- readRDS(here("processed-data", "02_SPEAQeasy", "rse_tx_living_brain_reanalysis_n516.Rds"))

## phenotype data
rse_tx_full$COI <- factor(ifelse(rse_tx_full$isPostMortem, "PM", "LIV"))
rse_tx_full$postmortem <- as.numeric(rse_tx_full$COI) - 1



## read in degradation data
rse_tx_degrade_full <- readRDS(file.path(dir_rdata, "rse_degrade_tx.Rds"))


table(rownames(rse_tx_full) %in% rownames(rse_tx_degrade_full))
# FALSE   TRUE
#  1655 196438

## Subset to the same transcripts measured here
rse_tx <- rse_tx_full[rownames(rse_tx_full) %in% rownames(rse_tx_degrade_full), ]
rse_tx_degrade <- rse_tx_degrade_full[rownames(rse_tx), ]

## Further subset the degradation data to just the RiboZero samples
rse_tx_degrade <- rse_tx_degrade[, rse_tx_degrade$sra_attribute.library_type == "RiboZero"]

## Also update the annotation information
rowRanges(rse_tx_degrade) <- rowRanges(rse_tx)

## Fix some variables in rse_tx_degrade
rse_tx_degrade$BrNum <- factor(rse_tx_degrade$sra_attribute.brain_number)
rse_tx_degrade$degrade_time <- as.numeric(rse_tx_degrade$sra_attribute.degradation_time)

## Sort by brain donor and degradation time
rse_tx_degrade <- rse_tx_degrade[, order(rse_tx_degrade$BrNum, rse_tx_degrade$degrade_time)]
as.data.frame(colData(rse_tx_degrade)[, c("BrNum", "degrade_time")])
#             BrNum degrade_time
# SRR5643537 Br1385            0
# SRR5643538 Br1385           15
# SRR5643535 Br1385           30
# SRR5643536 Br1385           60
# SRR5643533 Br1729            0
# SRR5643534 Br1729           15
# SRR5643531 Br1729           30
# SRR5643532 Br1729           60
# SRR5643539 Br2015            0
# SRR5643540 Br2015           15
# SRR5643548 Br2015           30
# SRR5643547 Br2015           60
# SRR5643530 Br2020            0
# SRR5643529 Br2020           15
# SRR5643544 Br2020           30
# SRR5643543 Br2020           60
# SRR5643546 Br2074            0
# SRR5643545 Br2074           15
# SRR5643542 Br2074           30
# SRR5643541 Br2074           60

## Extract the tpm values
tpm <- assay(rse_tx, "tpm")
tpm_degrade <- assay(rse_tx_degrade, "TPM")

### recalculate degradation stas for all isoforms, not just 45k in objects
fit_degrade <- lmFit(log2(tpm_degrade + 1), model.matrix(~ degrade_time + BrNum, data = colData(rse_tx_degrade)))
degradation_tstats <- topTable(
    eBayes(fit_degrade),
    coef = 2,
    sort = "none",
    n = nrow(tpm_degrade),
    genelist = rowData(rse_tx_degrade)
)
degradation_tstats <- degradation_tstats[order(degradation_tstats$P.Value), ]


## plots for top degradation t-stats
# filter to rows
cell_rows <- rownames(tpm) %in% transcripts$cell_component
sum(cell_rows)
# [1] 2924
standard_rows <- rownames(tpm) %in% transcripts$standard
sum(standard_rows)
# [1] 1760

top10_rows <- rownames(tpm) %in% rownames(degradation_tstats)[1:10]
top25_rows <- rownames(tpm) %in% rownames(degradation_tstats)[1:25]
top50_rows <- rownames(tpm) %in% rownames(degradation_tstats)[1:50]
top100_rows <- rownames(tpm) %in% rownames(degradation_tstats)[1:100]
top200_rows <- rownames(tpm) %in% rownames(degradation_tstats)[1:200]

## make list of TPMs
tpm_list <- list(
    tpm_standard = tpm[standard_rows, ],
    tpm_standard_pm = tpm[standard_rows, rse_tx$COI == "PM"],
    tpm_standard_liv = tpm[standard_rows, rse_tx$COI == "LIV"],
    tpm_cell = tpm[cell_rows, ],
    tpm_cell_pm = tpm[cell_rows, rse_tx$COI == "PM"],
    tpm_cell_liv = tpm[cell_rows, rse_tx$COI == "LIV"],
    tpm_top10 = tpm[top10_rows, ],
    tpm_top25 = tpm[top25_rows, ],
    tpm_top50 = tpm[top50_rows, ],
    tpm_top100 = tpm[top100_rows, ],
    tpm_top200 = tpm[top200_rows, ]
)
lapply(tpm_list, dim)
# $tpm_standard
# [1] 1760  516
#
# $tpm_standard_pm
# [1] 1760  243
#
# $tpm_standard_liv
# [1] 1760  273
#
# $tpm_cell
# [1] 2924  516
#
# $tpm_cell_pm
# [1] 2924  243
#
# $tpm_cell_liv
# [1] 2924  273
#
# $tpm_top10
# [1]  10 516
#
# $tpm_top25
# [1]  25 516
#
# $tpm_top50
# [1]  50 516
#
# $tpm_top100
# [1] 100 516
#
# $tpm_top200
# [1] 200 516

set.seed(743)
qsva_list <- mclapply(tpm_list, function(x) {
    expr <- log2(x + 1)
    m <- matrix(1, nrow = ncol(expr), ncol = 1)
    k <- sva::num.sv(expr, m)
    pca <- prcomp(t(expr))
    qsvs <- pca$x[, 1:k]
    var_expl <- getPcaVars(pca)[1:k]
    list(qSVs = qsvs, varExpl = var_expl)
}, mc.cores = 6)

## Save for later
save(qsva_list, file = file.path(dir_rdata, "qsva_objects.RData"))


#############
## load objects
# load(file.path(dir_rdata, "qsva_objects.RData"), verbose = TRUE)

## write out as text
for (i in seq(along = qsva_list)) {
    write.csv(qsva_list[[i]]$qSVs,
        file = file.path(dir_rdata, paste0(
            "qSVs_",
            gsub("tpm_", "", names(qsva_list)[i]), ".csv"
        ))
    )
}

# associations to PM vs LIV
coi_assoc <- lapply(qsva_list[c(1, 4, 7:11)], function(y) {
    round(-log10(apply(
        t(t(y$qSVs)), 2,
        function(x) {
            t.test(x ~ rse_tx$COI)$p.value
        }
    )), 1)
})
coi_assoc
# $tpm_standard
#  PC1  PC2  PC3  PC4  PC5  PC6  PC7  PC8  PC9 PC10
# 26.5  0.6 50.3  3.3 18.3  0.0  5.9  1.3  0.8  0.6
#
# $tpm_cell
#  PC1  PC2  PC3  PC4  PC5  PC6  PC7  PC8  PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16
# 35.0  8.0 26.2  2.1 26.9  1.1  0.5  0.9  1.6  0.0  0.3  0.5  0.5  6.7  0.3  0.8
#
# $tpm_top10
# [1] 8
#
# $tpm_top25
# [1] 8.6
#
# $tpm_top50
# [1] 18.6
#
# $tpm_top100
#  PC1  PC2  PC3
# 21.1 62.2  0.7
#
# $tpm_top200
#  PC1  PC2  PC3
# 21.0 13.8 44.8

degradation_tstats[1:10, ]
#                   source       type score phase            gene_id                          gene_type gene_status     gene_name level           havana_gene     transcript_id      transcript_type
# ENST00000623212.1 HAVANA transcript    NA    NA  ENSG00000278935.1                                TEC       KNOWN  RP11-173D3.4     2  OTTHUMG00000175828.1 ENST00000623212.1                  TEC
# ENST00000507754.8 HAVANA transcript    NA    NA ENSG00000186010.18                     protein_coding       KNOWN       NDUFA13     2  OTTHUMG00000162211.8 ENST00000507754.8       protein_coding
# ENST00000602936.1 HAVANA transcript    NA    NA  ENSG00000233012.2   transcribed_processed_pseudogene       KNOWN       HDAC1P2     2  OTTHUMG00000037360.2 ENST00000602936.1 processed_transcript
# ENST00000636670.1 HAVANA transcript    NA    NA ENSG00000141837.19                     protein_coding       KNOWN       CACNA1A     2 OTTHUMG00000044590.16 ENST00000636670.1      retained_intron
# ENST00000540997.2 HAVANA transcript    NA    NA ENSG00000185495.10 transcribed_unprocessed_pseudogene       KNOWN RP11-504P24.3     2  OTTHUMG00000037472.2 ENST00000540997.2      retained_intron
# ENST00000489327.1 HAVANA transcript    NA    NA ENSG00000115977.18                     protein_coding       KNOWN          AAK1     2  OTTHUMG00000129648.9 ENST00000489327.1      retained_intron
# ENST00000381759.8 HAVANA transcript    NA    NA ENSG00000127663.14                     protein_coding       KNOWN         KDM4B     2  OTTHUMG00000180281.4 ENST00000381759.8       protein_coding
# ENST00000429684.1 HAVANA transcript    NA    NA  ENSG00000225813.1               processed_pseudogene       KNOWN    AC009299.4     2  OTTHUMG00000153884.1 ENST00000429684.1 processed_pseudogene
# ENST00000457898.1 HAVANA transcript    NA    NA  ENSG00000226029.1                            lincRNA       KNOWN  RP4-798A10.2     2  OTTHUMG00000002314.1 ENST00000457898.1              lincRNA
# ENST00000301894.6 HAVANA transcript    NA    NA ENSG00000110076.18                     protein_coding       KNOWN         NRXN2     2 OTTHUMG00000045214.18 ENST00000301894.6       protein_coding
#                   transcript_status   transcript_name transcript_support_level                    tag    havana_transcript exon_number exon_id         ont        protein_id      ccdsid
# ENST00000623212.1             KNOWN  RP11-173D3.4-001                       NA                  basic OTTHUMT00000431140.1        <NA>    <NA>        <NA>              <NA>        <NA>
# ENST00000507754.8             KNOWN       NDUFA13-001                        1                   CCDS OTTHUMT00000367916.6        <NA>    <NA>        <NA> ENSP00000423673.1 CCDS12404.2
# ENST00000602936.1             KNOWN       HDAC1P2-002                       NA                  basic OTTHUMT00000467811.1        <NA>    <NA>        <NA>              <NA>        <NA>
# ENST00000636670.1             KNOWN       CACNA1A-056                     <NA> RNA_Seq_supported_only OTTHUMT00000489462.1        <NA>    <NA>        <NA>              <NA>        <NA>
# ENST00000540997.2             KNOWN RP11-504P24.3-006                       NA                  basic OTTHUMT00000478307.1        <NA>    <NA>        <NA>              <NA>        <NA>
# ENST00000489327.1             KNOWN          AAK1-005                        5                   <NA> OTTHUMT00000327346.1        <NA>    <NA>        <NA>              <NA>        <NA>
# ENST00000381759.8             KNOWN         KDM4B-011                        1                  basic OTTHUMT00000450559.1        <NA>    <NA>        <NA> ENSP00000371178.3        <NA>
# ENST00000429684.1             KNOWN    AC009299.4-001                       NA                  basic OTTHUMT00000332841.1        <NA>    <NA> PGO:0000004              <NA>        <NA>
# ENST00000457898.1             KNOWN  RP4-798A10.2-001                        2                  basic OTTHUMT00000006685.1        <NA>    <NA>        <NA>              <NA>        <NA>
# ENST00000301894.6             KNOWN         NRXN2-003                        5                   CCDS OTTHUMT00000141952.1        <NA>    <NA>        <NA> ENSP00000301894.2  CCDS8078.1
#                         logFC   AveExpr         t      P.Value    adj.P.Val        B
# ENST00000623212.1  0.01452433 0.5526962  16.03968 1.205789e-10 6.249316e-06 12.38480
# ENST00000507754.8 -0.01902535 5.0288771 -16.01634 1.230275e-10 6.249316e-06 12.36339
# ENST00000602936.1  0.01856818 1.0348904  15.79744 1.487400e-10 6.249316e-06 12.16128
# ENST00000636670.1  0.03684928 1.6021253  15.72812 1.580311e-10 6.249316e-06 12.09676
# ENST00000540997.2  0.01946838 1.2727543  15.40293 2.106675e-10 6.249316e-06 11.79064
# ENST00000489327.1  0.02478279 1.3632369  15.26295 2.388181e-10 6.249316e-06 11.65711
# ENST00000381759.8  0.01748177 0.7741020  15.07536 2.829856e-10 6.249316e-06 11.47646
# ENST00000429684.1  0.03381886 2.2386163  14.80423 3.628540e-10 6.249316e-06 11.21182
# ENST00000457898.1  0.01719272 1.2009255  14.77432 3.730345e-10 6.249316e-06 11.18237
# ENST00000301894.6  0.03929402 2.1438709  14.72692 3.897966e-10 6.249316e-06 11.13559

pdf(file.path(dir_plots, "qSV1_top10_vs_COI.pdf"))
boxplot(
    qsva_list$tpm_top10$qSVs ~ rse_tx$COI,
    ylab = paste0(
        "qSV1 (top10 Tx: ", qsva_list$tpm_top10$varExpl[1],
        "% Var Expl)"
    ),
    xlab = ""
)
dev.off()


### check qSV vs COI
apply(qsva_list$tpm_standard$qSVs, 2, function(x) {
    t.test(x ~ rse_tx$COI)$p.value
})
#          PC1          PC2          PC3          PC4          PC5          PC6          PC7          PC8          PC9         PC10
# 3.130702e-27 2.474994e-01 4.735134e-51 4.858927e-04 4.986386e-19 9.987815e-01 1.330729e-06 5.471337e-02 1.682187e-01 2.552738e-01
apply(qsva_list$tpm_cell$qSVs, 2, function(x) {
    t.test(x ~ rse_tx$COI)$p.value
})
#          PC1          PC2          PC3          PC4          PC5          PC6          PC7          PC8          PC9         PC10         PC11         PC12
# 1.049906e-35 9.237539e-09 6.699887e-27 7.276874e-03 1.313578e-27 8.246924e-02 3.104048e-01 1.292980e-01 2.785985e-02 9.449539e-01 4.622598e-01 2.970995e-01
#         PC13         PC14         PC15         PC16
# 3.416701e-01 2.116482e-07 5.622145e-01 1.510465e-01

## top qSV is very correlated
cor(qsva_list$tpm_standard$qSVs[, 1], qsva_list$tpm_cell$qSVs[, 1])
# [1] 0.9943452

qsva_list$tpm_standard$varExpl
# [1] 63.200  6.170  5.020  3.270  2.920  1.620  1.220  1.030  0.917  0.681

pdf(file.path(dir_plots, "qSV1_standard_vs_COI.pdf"))
boxplot(
    qsva_list$tpm_standard$qSVs[, 1] ~ rse_tx$COI,
    ylab = paste0(
        "qSV1 (standard: ",
        qsva_list$tpm_standard$varExpl[1],
        "% Var Expl)"
    ),
    xlab = ""
)
dev.off()

pdf(file.path(dir_plots, "qSV1_cell_vs_COI.pdf"))
boxplot(
    qsva_list$tpm_cell$qSVs[, 1] ~ rse_tx$COI,
    ylab = paste0(
        "qSV1 (cell: ",
        qsva_list$tpm_cell$varExpl[1],
        "% Var Expl)"
    ),
    xlab = ""
)
dev.off()


## Compare against RNAseq metrics
gIndexes <- splitit(rse_tx$COI)

pdf(file.path(dir_plots, "qSV1_top10_vs_totalAssignedGene.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 2,
    cex.lab = 2
)
palette(brewer.pal(4, "Dark2"))
plot(
    qsva_list$tpm_top10$qSVs,
    rse_tx$totalAssignedGene,
    xlab = paste0(
        "qSV1 (top10 Tx: ", qsva_list$tpm_top10$varExpl[1],
        "% Var Expl)"
    ),
    ylab = "totalAssignedGene",
    pch = 21,
    bg = rse_tx$COI
)
for (i in seq(along = gIndexes)) {
    ii <- gIndexes[[i]]
    abline(lm(rse_tx$totalAssignedGene ~ qsva_list$tpm_top10$qSVs, subset = ii),
        lwd = 4,
        col = i
    )
}
dev.off()


pdf(file.path(dir_plots, "qSV1_top10_vs_mitoRate.pdf"))
par(
    mar = c(4, 6, 2, 2),
    cex.axis = 2,
    cex.lab = 2
)
palette(brewer.pal(4, "Dark2"))
plot(
    qsva_list$tpm_top10$qSVs,
    rse_tx$mitoRate,
    xlab = paste0(
        "qSV1 (top10 Tx: ", qsva_list$tpm_top10$varExpl[1],
        "% Var Expl)"
    ),
    ylab = "mitoRate",
    pch = 21,
    bg = rse_tx$COI
)
for (i in seq(along = gIndexes)) {
    ii <- gIndexes[[i]]
    abline(lm(rse_tx$mitoRate ~ qsva_list$tpm_top10$qSVs, subset = ii),
        lwd = 4,
        col = i
    )
}
dev.off()

###########
## de checks

## For simplicity later on, I'll use this variable name
rse_tx$samplingAge <- rse_tx$ageDeath_num.biospecimen

## Fix race since it has 206 NAs
rse_tx$race <- as.character(rse_tx$race)
rse_tx$race[is.na(rse_tx$race)] <- "Unknown"

## add cell comp PCs
cellProps <- read.csv(here("processed-data", "05_check_gene_exprs", "LBP_burkeDecon.csv"), row.names = 1)
cellProps <- cellProps[colnames(rse_gene), ]
cellPca <- prcomp(cellProps)
getPcaVars(cellPca)[1:2]
# [1] 79.0 11.1
rse_tx$cellPC <- cellPca$x[, 1]

## model
mod <- model.matrix(
    ~ COI + race + sex + diagnosis +
        RIN + cellPC +
        mitoRate +
        totalAssignedGene +
        rRNA_rate +
        overallMapRate,
    data = colData(rse_tx)
)

## modeling
f <- lmFit(log2(tpm + 1), mod)
out <- topTable(
    eBayes(f),
    coef = 2,
    sort = "none",
    n = nrow(tpm),
    genelist = rowData(rse_tx)
)

## against degradation
out_top10000 <- out[order(out$P.Value)[1:10000], ]
degradation_tstats_top10000 <- degradation_tstats[match(
    out_top10000$transcript_id,
    degradation_tstats$transcript_id
), ]
pdf(file.path(dir_plots, "top10k_t_vs_degradation_t.pdf"))
plot(out_top10000$t, degradation_tstats_top10000$t)
dev.off()
tt_10000 <- table(`PM>LIV` = out_top10000$t > 0, `Degrade>0` = degradation_tstats_top10000$t > 0)
tt_10000
#        Degrade>0
# PM>LIV  FALSE TRUE
#   FALSE    11   14
#   TRUE    920 9055
getOR(tt_10000)
# [1] 7.733307

#####
f_standard <- lmFit(log2(tpm + 1), cbind(mod, qsva_list$tpm_standard$qSVs))
out_standard <- topTable(
    eBayes(f_standard),
    coef = 2,
    sort = "none",
    n = nrow(tpm),
    genelist = rowData(rse_tx)
)

f_cell <- lmFit(log2(tpm + 1), cbind(mod, qsva_list$tpm_cell$qSVs))
out_cell <- topTable(
    eBayes(f_cell),
    coef = 2,
    sort = "none",
    n = nrow(tpm),
    genelist = rowData(rse_tx)
)

cor(out_cell$t, out$t)
# [1] 0.1251368
pdf(file.path(dir_plots, "smootScatter_out_cell_t_vs_out_t.pdf"))
smoothScatter(out_cell$t, out$t)
dev.off()

cor(out_standard$t, out$t)
# [1] 0.09316496
pdf(file.path(dir_plots, "smootScatter_out_standard_t_vs_out_t.pdf"))
smoothScatter(out_standard$t, out$t)
dev.off()

cor(out_cell$t, out_standard$t)
# [1] 0.6607493
pdf(file.path(dir_plots, "smootScatter_out_cell_t_vs_out_standard_t.pdf"))
smoothScatter(out_cell$t, out_standard$t)
dev.off()

out_cell[order(out_cell$P.Value)[1:25], c("gene_name", "t", "P.Value")]
#                       gene_name          t      P.Value
# ENST00000563103.1 CTD-2651B20.6 -17.244021 3.012677e-52
# ENST00000585776.5  RP11-383M4.6  12.950420 4.001028e-33
# ENST00000385129.3        MIR27B -12.623045 9.281923e-32
# ENST00000411199.1     RNA5SP304 -12.343306 1.320583e-30
# ENST00000517961.3    AC084082.3 -11.716965 4.508522e-28
# ENST00000546283.5        NDUFS7 -10.264858 1.727036e-22
# ENST00000305943.7        RNF187 -10.109974 6.403418e-22
# ENST00000412681.2          NRGN -10.100416 6.939923e-22
# ENST00000601154.1         CCDC9   9.944578 2.559176e-21
# ENST00000527676.2          ETS1   9.556798 6.216946e-20
# ENST00000546428.2      SLC25A21   9.457360 1.389995e-19
# ENST00000568314.1 CTD-2651B20.7  -9.411831 2.005386e-19
# ENST00000362263.3      MIRLET7D  -9.385210 2.483319e-19
# ENST00000394067.6          KLC2  -9.324177 4.047737e-19
# ENST00000558365.1     ADAMTS7P4   9.316018 4.320203e-19
# ENST00000506696.1          SNCB  -9.223221 9.039253e-19
# ENST00000367906.7          RGS4  -9.043146 3.733165e-18
# ENST00000239223.3         DUSP1   8.958207 7.239665e-18
# ENST00000310112.7          SNCB  -8.937419 8.508062e-18
# ENST00000590926.1         APLP1  -8.922216 9.572649e-18
# ENST00000315707.3     LINC00324  -8.752419 3.538215e-17
# ENST00000335312.7       PIP5K1C  -8.749487 3.618438e-17
# ENST00000624231.1 RP11-114G11.4   8.631707 8.866434e-17
# ENST00000414364.1     KIF25-AS1   8.571490 1.397360e-16
# ENST00000553375.1       ZFP36L1   8.309787 9.826687e-16

f_top100 <- lmFit(log2(tpm + 1), cbind(mod, qsva_list$tpm_top100$qSVs))
out_top100 <- topTable(
    eBayes(f_top100),
    coef = 2,
    sort = "none",
    n = nrow(tpm),
    genelist = rowData(rse_tx)
)

table(out$adj.P.Val < 0.05)
# FALSE   TRUE
# 84546 111892
table(out_standard$adj.P.Val < 0.05)
#  FALSE   TRUE
# 175870  20568
table(out_cell$adj.P.Val < 0.05)
#  FALSE   TRUE
# 191289   5149
table(out_top100$adj.P.Val < 0.05)
#  FALSE   TRUE
# 154976  41462

degradation_tstats_out <- degradation_tstats[match(out$transcript_id, degradation_tstats$transcript_id), ]

pdf(file.path(dir_plots, "smootScatter_out_t_vs_degradation_t.pdf"))
smoothScatter(out$t, degradation_tstats_out$t)
dev.off()

sumt <- degradation_tstats_out$t + out$t
# cir_ind = c(which.min(sumt), which.max(sumt))

pdf(file.path(dir_plots, "Dequal_txLevel.pdf"))
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
plot(
    out$t,
    degradation_tstats_out$t,
    pch = 21,
    bg = "grey",
    xlab = "PM vs LIV (t-stat)",
    ylab = "Degradation (t-stat)"
)
# points(out$t[cir_ind], degradation_tstats_out$t[cir_ind],
# 	pch = 1, col = c("red","blue"), cex=1.5)
legend("topleft", paste0(
    "Cor=",
    signif(
        cor(out$t, degradation_tstats_out$t, use = "comp"), 3
    )
),
cex = 1.5
)
abline(h = 0, v = 0, lty = 2)
dev.off()


# ######################
#### example plots #####
plot_ind <- c(order(sumt)[1:25], order(sumt, decreasing = TRUE)[1:25])
degradation_tstats_sumt <- degradation_tstats_out[plot_ind, ]
tpm_degrade_order <- as.matrix(tpm_degrade[rownames(degradation_tstats_sumt), ])
bIndexes <- splitit(rse_tx_degrade$BrNum)

pdf(file.path(dir_plots, "degradation_transcripts_top.pdf"))
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
palette(brewer.pal(5, "Set1"))
for (i in 1:nrow(tpm_degrade_order)) {
    plot(
        log2(tpm_degrade_order[i, ] + 1) ~ rse_tx_degrade$degrade_time,
        pch = 21,
        bg = rse_tx_degrade$BrNum,
        cex = 2,
        xlab = "Degradation Time (min)",
        ylab = "log2(TPM+1)",
        main = paste0(
            rownames(degradation_tstats_sumt)[i],
            "\n",
            degradation_tstats_sumt$gene_name[i],
            " (",
            degradation_tstats_sumt$transcript_type[i],
            ")"
        )
    )
    for (j in seq(along = bIndexes)) {
        lines(
            log2(tpm_degrade_order[i, ] + 1) ~ rse_tx_degrade$degrade_time,
            lwd = 2,
            col = j,
            subset = bIndexes[[j]]
        )
    }
    lc <- ifelse(degradation_tstats_sumt$t[i] > 0, "topleft", "topright")
    legend(lc, paste0("p=", signif(degradation_tstats_sumt$P.Value[i], 3)),
        cex =
            1.7
    )
}
dev.off()

tpm_lbp_order <- as.matrix(tpm[rownames(tpm_degrade_order), ])
pdf(file.path(dir_plots, "degradation_transcripts_lbp.pdf"))
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
palette(brewer.pal(5, "Dark2"))
for (i in 1:nrow(tpm_degrade_order)) {
    boxplot(
        log2(tpm_lbp_order[i, ] + 1) ~ rse_tx$COI,
        xlab = "",
        ylab = "log2(TPM+1)",
        ylim = range(log2(tpm_lbp_order[i, ] + 1)),
        outline = FALSE,
        main = paste0(
            rownames(degradation_tstats_sumt)[i],
            "\n",
            degradation_tstats_sumt$gene_name[i],
            " (",
            degradation_tstats_sumt$transcript_type[i],
            ")"
        )
    )

    points(
        log2(tpm_lbp_order[i, ] + 1) ~ jitter(as.numeric(rse_tx$COI), amount = 0.15),
        pch = 21,
        bg = rse_tx$COI
    )
}
dev.off()

### prediction modeling
X <- as.data.frame(t(tpm_degrade[rownames(degradation_tstats)[1:10], ]))
X_new <- as.data.frame(t(tpm[colnames(X), ]))
X$Time <- rse_tx_degrade$degrade_time

f_pred <- lm(Time ~ ., data = X)

pred_time <- predict(f_pred, X_new)

pdf(file.path(dir_plots, "degradation_prediction_lbp.pdf"))
par(
    mar = c(5, 6, 2, 2),
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
palette(brewer.pal(5, "Dark2"))
boxplot(
    pred_time ~ rse_tx$COI,
    xlab = "",
    ylab = "Predicted Degradation Time (Min)",
    ylim = range(pred_time),
    outline = FALSE
)
points(pred_time ~ jitter(as.numeric(rse_tx$COI), amount = 0.15),
    pch = 21,
    bg = rse_tx$COI
)
legend("topleft", paste0(
    "p=",
    signif(summary(
        lm(pred_time ~ rse_tx$COI)
    )$coef[2, 4], 3)
), cex = 2)
dev.off()

##########
## age + sex + dx?
##########

## in living
pheno_liv <- colData(rse_tx)[rse_tx$COI == "LIV", ]
pheno_liv$diagnosis <- droplevels(pheno_liv$diagnosis)
mod_liv <- model.matrix(~ samplingAge + sex + diagnosis + race + rnaBatch + libraryBatch,
    data = pheno_liv
)
tpm_liv <- tpm[, rownames(mod_liv)]
standard_qsv_liv <- qsva_list$tpm_standard_liv$qSVs[rownames(mod_liv), ]
cell_qsv_liv <- qsva_list$tpm_cell_liv$qSVs[rownames(mod_liv), ]

f_liv <- lmFit(log2(tpm_liv + 1), mod_liv)
f_liv_standard <- lmFit(log2(tpm_liv + 1), cbind(mod_liv, standard_qsv_liv))
f_liv_cell <- lmFit(log2(tpm_liv + 1), cbind(mod_liv, cell_qsv_liv))

## in pm
pheno_pm <- colData(rse_tx)[rse_tx$COI == "PM", ]
pheno_pm$diagnosis <- droplevels(pheno_pm$diagnosis)
mod_pm <- model.matrix(~ samplingAge + sex + diagnosis + race + rnaBatch + libraryBatch,
    data = pheno_pm
)
tpm_pm <- tpm[, rownames(mod_pm)]
standard_qsv_pm <- qsva_list$tpm_standard_pm$qSVs[rownames(mod_pm), ]
cell_qsv_pm <- qsva_list$tpm_cell_pm$qSVs[rownames(mod_pm), ]

f_pm <- lmFit(log2(tpm_pm + 1), mod_pm)
f_pm_standard <- lmFit(log2(tpm_pm + 1), cbind(mod_pm, standard_qsv_pm))
f_pm_cell <- lmFit(log2(tpm_pm + 1), cbind(mod_pm, cell_qsv_pm))

#############
##### age ###
age_liv <- topTable(eBayes(f_liv),
    coef = 2,
    sort = "none",
    n = nrow(tpm)
)
age_liv_standard <- topTable(
    eBayes(f_liv_standard),
    coef = 2,
    sort = "none",
    n = nrow(tpm)
)
age_liv_cell <- topTable(eBayes(f_liv_cell),
    coef = 2,
    sort = "none",
    n = nrow(tpm)
)
age_pm <- topTable(eBayes(f_pm),
    coef = 2,
    sort = "none",
    n = nrow(tpm)
)
age_pm_standard <- topTable(
    eBayes(f_pm_standard),
    coef = 2,
    sort = "none",
    n = nrow(tpm)
)
age_pm_cell <- topTable(eBayes(f_pm_cell),
    coef = 2,
    sort = "none",
    n = nrow(tpm)
)

## liv discovery
mean(age_pm$P.Value[age_liv$P.Value < 1e-4] < 0.05 &
    sign(age_pm$t[age_liv$P.Value < 1e-4]) ==
        sign(age_liv$t[age_liv$P.Value < 1e-4]))
# [1] 0.2121212
mean(
    age_pm_standard$P.Value[age_liv_standard$P.Value < 1e-4] < 0.05 &
        sign(age_pm_standard$t[age_liv_standard$P.Value < 1e-4]) ==
            sign(age_liv_standard$t[age_liv_standard$P.Value < 1e-4])
)
# [1] 0.257764
mean(age_pm_cell$P.Value[age_liv_cell$P.Value < 1e-4] < 0.05 &
    sign(age_pm_cell$t[age_liv_cell$P.Value < 1e-4]) ==
        sign(age_liv_cell$t[age_liv_cell$P.Value < 1e-4]))
# [1] 0.2765957

## pm discovery
mean(age_liv$P.Value[age_pm$P.Value < 1e-4] < 0.05 &
    sign(age_liv$t[age_pm$P.Value < 1e-4]) ==
        sign(age_pm$t[age_pm$P.Value < 1e-4]))
# [1] 0.0115138
mean(
    age_liv_standard$P.Value[age_pm_standard$P.Value < 1e-4] < 0.05 &
        sign(age_liv_standard$t[age_pm_standard$P.Value < 1e-4]) ==
            sign(age_pm_standard$t[age_pm_standard$P.Value < 1e-4])
)
# [1] 0.2172619
mean(age_liv_cell$P.Value[age_pm_cell$P.Value < 1e-4] < 0.05 &
    sign(age_liv_cell$t[age_pm_cell$P.Value < 1e-4]) ==
        sign(age_pm_cell$t[age_pm_cell$P.Value < 1e-4]))
# [1] 0.2823529

pdf(file.path(dir_plots, "smootScatter_age_liv_t_vs_age_pm_t.pdf"))
smoothScatter(age_liv$t, age_pm$t)
dev.off()
pdf(file.path(dir_plots, "smootScatter_age_liv_standard_t_vs_age_pm_standard_t.pdf"))
smoothScatter(age_liv_standard$t, age_pm_standard$t)
dev.off()
pdf(file.path(dir_plots, "smootScatter_age_liv_cell_t_vs_age_pm_cell_t.pdf"))
smoothScatter(age_liv_cell$t, age_pm_cell$t)
dev.off()

#############
##### sex ###
colnames(coef(eBayes(f_liv)))
#  [1] "(Intercept)"                            "samplingAge"                            "sexmale"                                "diagnosisessential tremor"
#  [5] "diagnosisMajor Depressive Disorder"     "diagnosisobsessive compulsive disorder" "diagnosisParkinson's disease"           "raceBlack or African American"
#  [9] "raceUnknown"                            "raceWhite"                              "rnaBatchLBPSEMA4_PLATE_2"               "rnaBatchLBPSEMA4_PLATE_3"
# [13] "rnaBatchLBPSEMA4_PLATE_4"               "rnaBatchLBPSEMA4_PLATE_5"               "rnaBatchLBPSEMA4_PLATE_6"               "rnaBatchLBPSEMA4_PLATE_7"
# [17] "rnaBatchLBPSEMA4_PLATE_8"               "rnaBatchLBPSEMA4_PLATE_9"               "libraryBatchdb10"                       "libraryBatchdb11"
# [21] "libraryBatchdb12"                       "libraryBatchdb13"                       "libraryBatchdb14"                       "libraryBatchdb15"
# [25] "libraryBatchdb16"                       "libraryBatchdb17"                       "libraryBatchdb18"                       "libraryBatchdb19"
# [29] "libraryBatchdb2"                        "libraryBatchdb20"                       "libraryBatchdb21"                       "libraryBatchdb22"
# [33] "libraryBatchdb23"                       "libraryBatchdb24"                       "libraryBatchdb25"                       "libraryBatchdb3"
# [37] "libraryBatchdb4"                        "libraryBatchdb5"                        "libraryBatchdb6"                        "libraryBatchdb7"
# [41] "libraryBatchdb8"                        "libraryBatchdb9"
sex_liv <- topTable(eBayes(f_liv),
    coef = "sexmale",
    sort = "none",
    n = nrow(tpm)
)
sex_liv_standard <- topTable(
    eBayes(f_liv_standard),
    coef = "sexmale",
    sort = "none",
    n = nrow(tpm)
)
sex_liv_cell <- topTable(eBayes(f_liv_cell),
    coef = "sexmale",
    sort = "none",
    n = nrow(tpm)
)
sex_pm <- topTable(eBayes(f_pm),
    coef = "sexmale",
    sort = "none",
    n = nrow(tpm)
)
sex_pm_standard <- topTable(
    eBayes(f_pm_standard),
    coef = "sexmale",
    sort = "none",
    n = nrow(tpm)
)
sex_pm_cell <- topTable(eBayes(f_pm_cell),
    coef = "sexmale",
    sort = "none",
    n = nrow(tpm)
)

## liv discovery
mean(sex_pm$P.Value[sex_liv$P.Value < 1e-4] < 0.05 &
    sign(sex_pm$t[sex_liv$P.Value < 1e-4]) ==
        sign(sex_liv$t[sex_liv$P.Value < 1e-4]))
# [1] 0.941896
mean(
    sex_pm_standard$P.Value[sex_liv_standard$P.Value < 1e-4] < 0.05 &
        sign(sex_pm_standard$t[sex_liv_standard$P.Value < 1e-4]) ==
            sign(sex_liv_standard$t[sex_liv_standard$P.Value < 1e-4])
)
# [1] 0.857732
mean(sex_pm_cell$P.Value[sex_liv_cell$P.Value < 1e-4] < 0.05 &
    sign(sex_pm_cell$t[sex_liv_cell$P.Value < 1e-4]) ==
        sign(sex_liv_cell$t[sex_liv_cell$P.Value < 1e-4]))
# [1] 0.1496063

## pm discovery
mean(sex_liv$P.Value[sex_pm$P.Value < 1e-4] < 0.05 &
    sign(sex_liv$t[sex_pm$P.Value < 1e-4]) ==
        sign(sex_pm$t[sex_pm$P.Value < 1e-4]))
# [1] 0.6862302
mean(
    sex_liv_standard$P.Value[sex_pm_standard$P.Value < 1e-4] < 0.05 &
        sign(sex_liv_standard$t[sex_pm_standard$P.Value < 1e-4]) ==
            sign(sex_pm_standard$t[sex_pm_standard$P.Value < 1e-4])
)
# [1] 0.9254079
mean(sex_liv_cell$P.Value[sex_pm_cell$P.Value < 1e-4] < 0.05 &
    sign(sex_liv_cell$t[sex_pm_cell$P.Value < 1e-4]) ==
        sign(sex_pm_cell$t[sex_pm_cell$P.Value < 1e-4]))
# [1] 0.07801418

pdf(file.path(dir_plots, "smootScatter_sex_liv_t_vs_sex_pm_t.pdf"))
smoothScatter(sex_liv$t, sex_pm$t)
dev.off()
pdf(file.path(dir_plots, "smootScatter_sex_liv_standard_t_vs_sex_pm_standard_t.pdf"))
smoothScatter(sex_liv_standard$t, sex_pm_standard$t)
dev.off()
pdf(file.path(dir_plots, "smootScatter_sex_liv_cell_t_vs_sex_pm_cell_t.pdf"))
smoothScatter(sex_liv_cell$t, sex_pm_cell$t)
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
#  date     2023-09-08
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
#  biocthis               1.11.4    2023-07-04 [1] Github (lcolladotor/biocthis@d062a71)
#  biomaRt                2.56.1    2023-06-11 [1] Bioconductor
#  Biostrings             2.68.1    2023-05-16 [1] Bioconductor
#  bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
#  bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
#  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  blob                   1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
#  brio                   1.1.3     2021-11-30 [1] CRAN (R 4.3.0)
#  bumphunter             1.42.0    2023-04-25 [1] Bioconductor
#  cachem                 1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  callr                  3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.1)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
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
#  edgeR                  3.42.4    2023-05-31 [1] Bioconductor
#  ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
#  fansi                  1.0.4     2023-01-22 [1] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
#  filelock               1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
#  foreach                1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
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
#  ggplot2                3.4.3     2023-08-14 [1] CRAN (R 4.3.0)
#  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
#  googledrive            2.1.1     2023-06-11 [1] CRAN (R 4.3.0)
#  gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
#  HDF5Array              1.28.1    2023-05-01 [1] Bioconductor
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
#  hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools              0.5.6     2023-08-10 [1] CRAN (R 4.3.0)
#  htmlwidgets            1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
#  httpuv                 1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
#  httr                   1.4.7     2023-08-15 [1] CRAN (R 4.3.0)
#  illuminaio             0.42.0    2023-05-08 [1] Bioconductor
#  IRanges              * 2.34.1    2023-07-02 [1] Bioconductor
#  iterators              1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
#  jaffelab             * 0.99.32   2023-05-06 [1] Github (LieberInstitute/jaffelab@7b7afe3)
#  KEGGREST               1.40.0    2023-04-25 [1] Bioconductor
#  KernSmooth             2.23-22   2023-07-10 [1] CRAN (R 4.3.0)
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
#  mclust                 6.0.0     2022-10-31 [1] CRAN (R 4.3.0)
#  memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
#  mgcv                   1.9-0     2023-07-11 [1] CRAN (R 4.3.0)
#  mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  minfi                  1.46.0    2023-05-08 [1] Bioconductor
#  miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
#  multtest               2.56.0    2023-05-08 [1] Bioconductor
#  munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
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
#  qsvaR                * 1.4.0     2023-05-08 [1] Bioconductor
#  quadprog               1.5-8     2019-11-20 [1] CRAN (R 4.3.0)
#  R.cache                0.16.0    2022-07-21 [1] CRAN (R 4.3.0)
#  R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
#  R.oo                   1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
#  R.utils                2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
#  rappdirs               0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
#  RColorBrewer         * 1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
#  RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
#  readr                  2.1.4     2023-02-10 [1] CRAN (R 4.3.0)
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
#  rtracklayer          * 1.60.1    2023-08-20 [1] Bioconductor
#  S4Arrays               1.0.5     2023-07-30 [1] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [1] Bioconductor
#  scales                 1.2.1     2022-08-20 [1] CRAN (R 4.3.0)
#  scrime                 1.3.5     2018-12-01 [1] CRAN (R 4.3.0)
#  segmented              1.6-4     2023-04-13 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
#  shiny                  1.7.5     2023-08-12 [1] CRAN (R 4.3.0)
#  siggenes               1.74.0    2023-05-08 [1] Bioconductor
#  sparseMatrixStats      1.12.2    2023-07-02 [1] Bioconductor
#  stringi                1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
#  stringr                1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
#  styler                 1.10.2    2023-08-29 [1] CRAN (R 4.3.0)
#  SummarizedExperiment * 1.30.2    2023-06-11 [1] Bioconductor
#  suncalc                0.5.1     2022-09-29 [1] CRAN (R 4.3.0)
#  survival               3.5-7     2023-08-14 [1] CRAN (R 4.3.0)
#  sva                    3.48.0    2023-04-25 [1] Bioconductor
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
#  XVector                0.40.0    2023-04-25 [1] Bioconductor
#  yaml                   2.3.7     2023-01-23 [1] CRAN (R 4.3.0)
#  zlibbioc               1.46.0    2023-04-25 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
