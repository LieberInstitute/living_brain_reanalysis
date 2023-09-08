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

boxplot(
    qsva_list$tpm_standard$qSVs[, 1] ~ rse_tx$COI,
    ylab = paste0(
        "qSV1 (standard: ",
        qsva_list$tpm_standard$varExpl[1],
        "% Var Expl)"
    ),
    xlab = ""
)
boxplot(qsva_list$tpm_cell[, 1] ~ rse_tx$COI,
    ylab = "qSV1 (cell)",
    xlab = ""
)


### read in rnaseq metrics
metrics <- read.csv("phenotype/LBP_merged_seqmetrics.csv")
metrics <- metrics[match(rse_tx$specimenID, metrics$external_id), ]
plot(
    qsva_list$tpm_top10$qSVs,
    metrics$recount_qc.intron_sum_perc
)
plot(
    qsva_list$tpm_top10$qSVs,
    metrics$recount_qc.exon_fc.unique_perc
)
plot(
    qsva_list$tpm_top10$qSVs,
    metrics$"recount_qc.bc_auc.all_reads_annotated_bases"
)

###########
## de checks
pheno$race[pheno$race == "White "] <- "White"
pheno$race[pheno$race == ""] <- "Unknown"
pheno$diagnosis[pheno$diagnosis == "Parkinson's disease "] <- "Parkinson's disease"
pheno$samplingAge[pheno$samplingAge == "89+"] <- 89
pheno$samplingAge <- as.numeric(pheno$samplingAge)
rownames(pheno) <- pheno$specimenID
pheno$diagnosis <- factor(pheno$diagnosis)
pheno$diagnosis <- relevel(pheno$diagnosis, "control")

pheno <- cbind(pheno, metrics)

## model
mod <- model.matrix(
    ~ COI + race + sex + diagnosis +
        rin + cellPC + recount_qc.star.uniquely_mapped_reads_perc +
        recount_qc.aligned_readsperc.chrm +
        recount_qc.aligned_readsperc.chrx +
        recount_qc.bc_auc.all_perc +
        recount_qc.star.perc_of_reads_mapped_to_multiple_loci +
        recount_qc.gene_fc.all_perc +
        recount_qc.junction_avg_coverage,
    data = pheno
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
    ss(out_top10000$transcript_id, "\\."),
    ss(rownames(degradation_tstats), "\\.")
), ]
plot(out_top10000$t, degradation_tstats_top10000$t)
tt_10000 <- table(`PM>LIV` = out_top10000$t > 0, `Degrade>0` = degradation_tstats_top10000$t > 0)
tt_10000
getOR(tt_10000) # 13.6


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
smoothScatter(out_cell$t, out$t)
out_cell[order(out_cell$P.Value)[1:25], c("gene_name", "t", "P.Value")]

f_top100 <- lmFit(log2(tpm + 1), cbind(mod, qsva_list$tpm_top100$qSVs))
out_top100 <- topTable(
    eBayes(f_top100),
    coef = 2,
    sort = "none",
    n = nrow(tpm),
    genelist = rowData(rse_tx)
)

table(out$adj.P.Val < 0.05)
table(out_standard$adj.P.Val < 0.05)
table(out_cell$adj.P.Val < 0.05)
table(out_top100$adj.P.Val < 0.05)

degradation_tstats_out <- degradation_tstats[match(out$transcript_id, rownames(degradation_tstats)), ]

smoothScatter(out$t, degradation_tstats_out$t)

sumt <- degradation_tstats_out$t + out$t
# cir_ind = c(which.min(sumt), which.max(sumt))

pdf("plots/Dequal_txLevel.pdf")
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
bIndexes <- splitit(pheno_degrade$BrNum)

pdf("plots/degradation_transcripts_top.pdf")
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
palette(brewer.pal(5, "Set1"))
for (i in 1:nrow(tpm_degrade_order)) {
    plot(
        log2(tpm_degrade_order[i, ] + 1) ~ pheno_degrade$degrade_time,
        pch = 21,
        bg = pheno_degrade$BrNum,
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
            log2(tpm_degrade_order[i, ] + 1) ~ pheno_degrade$degrade_time,
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
pdf("plots/degradation_transcripts_lbp.pdf")
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
palette(brewer.pal(5, "Dark2"))
for (i in 1:nrow(tpm_degrade_order)) {
    boxplot(
        log2(tpm_lbp_order[i, ] + 1) ~ pheno$COI,
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
        log2(tpm_lbp_order[i, ] + 1) ~ jitter(as.numeric(pheno$COI), amount = 0.15),
        pch = 21,
        bg = pheno$COI
    )
}
dev.off()

### prediction modeling
X <- as.data.frame(t(tpm_degrade[rownames(degradation_tstats)[1:10], ]))
X_new <- as.data.frame(t(tpm[colnames(X), ]))
X$Time <- pheno_degrade$degrade_time

f_pred <- lm(Time ~ ., data = X)

pred_time <- predict(f_pred, X_new)

pdf("plots/degradation_prediction_lbp.pdf")
par(
    mar = c(5, 6, 2, 2),
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
palette(brewer.pal(5, "Dark2"))
boxplot(
    pred_time ~ pheno$COI,
    xlab = "",
    ylab = "Predicted Degradation Time (Min)",
    ylim = range(pred_time),
    outline = FALSE
)
points(pred_time ~ jitter(as.numeric(pheno$COI), amount = 0.15),
    pch = 21,
    bg = pheno$COI
)
legend("topleft", paste0(
    "p=",
    signif(summary(
        lm(pred_time ~ pheno$COI)
    )$coef[2, 4], 3)
), cex = 2)
dev.off()

##########
## age + sex + dx?
##########

## in living
pheno_liv <- pheno[pheno$COI == "LIV", ]
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
pheno_pm <- pheno[pheno$COI == "PM", ]
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
mean(
    age_pm_standard$P.Value[age_liv_standard$P.Value < 1e-4] < 0.05 &
        sign(age_pm_standard$t[age_liv_standard$P.Value < 1e-4]) ==
            sign(age_liv_standard$t[age_liv_standard$P.Value < 1e-4])
)
mean(age_pm_cell$P.Value[age_liv_cell$P.Value < 1e-4] < 0.05 &
    sign(age_pm_cell$t[age_liv_cell$P.Value < 1e-4]) ==
        sign(age_liv_cell$t[age_liv_cell$P.Value < 1e-4]))

## pm discovery
mean(age_liv$P.Value[age_pm$P.Value < 1e-4] < 0.05 &
    sign(age_liv$t[age_pm$P.Value < 1e-4]) ==
        sign(age_pm$t[age_pm$P.Value < 1e-4]))
mean(
    age_liv_standard$P.Value[age_pm_standard$P.Value < 1e-4] < 0.05 &
        sign(age_liv_standard$t[age_pm_standard$P.Value < 1e-4]) ==
            sign(age_pm_standard$t[age_pm_standard$P.Value < 1e-4])
)
mean(age_liv_cell$P.Value[age_pm_cell$P.Value < 1e-4] < 0.05 &
    sign(age_liv_cell$t[age_pm_cell$P.Value < 1e-4]) ==
        sign(age_pm_cell$t[age_pm_cell$P.Value < 1e-4]))


smoothScatter(age_liv$t, age_pm$t)
smoothScatter(age_liv_standard$t, age_pm_standard$t)
smoothScatter(age_liv_cell$t, age_pm_cell$t)

#############
##### sex ###
sex_liv <- topTable(eBayes(f_liv),
    coef = 9,
    sort = "none",
    n = nrow(tpm)
)
sex_liv_standard <- topTable(
    eBayes(f_liv_standard),
    coef = 9,
    sort = "none",
    n = nrow(tpm)
)
sex_liv_cell <- topTable(eBayes(f_liv_cell),
    coef = 9,
    sort = "none",
    n = nrow(tpm)
)
sex_pm <- topTable(eBayes(f_pm),
    coef = 9,
    sort = "none",
    n = nrow(tpm)
)
sex_pm_standard <- topTable(
    eBayes(f_pm_standard),
    coef = 9,
    sort = "none",
    n = nrow(tpm)
)
sex_pm_cell <- topTable(eBayes(f_pm_cell),
    coef = 9,
    sort = "none",
    n = nrow(tpm)
)

## liv discovery
mean(sex_pm$P.Value[sex_liv$P.Value < 1e-4] < 0.05 &
    sign(sex_pm$t[sex_liv$P.Value < 1e-4]) ==
        sign(sex_liv$t[sex_liv$P.Value < 1e-4]))
mean(
    sex_pm_standard$P.Value[sex_liv_standard$P.Value < 1e-4] < 0.05 &
        sign(sex_pm_standard$t[sex_liv_standard$P.Value < 1e-4]) ==
            sign(sex_liv_standard$t[sex_liv_standard$P.Value < 1e-4])
)
mean(sex_pm_cell$P.Value[sex_liv_cell$P.Value < 1e-4] < 0.05 &
    sign(sex_pm_cell$t[sex_liv_cell$P.Value < 1e-4]) ==
        sign(sex_liv_cell$t[sex_liv_cell$P.Value < 1e-4]))

## pm discovery
mean(sex_liv$P.Value[sex_pm$P.Value < 1e-4] < 0.05 &
    sign(sex_liv$t[sex_pm$P.Value < 1e-4]) ==
        sign(sex_pm$t[sex_pm$P.Value < 1e-4]))
mean(
    sex_liv_standard$P.Value[sex_pm_standard$P.Value < 1e-4] < 0.05 &
        sign(sex_liv_standard$t[sex_pm_standard$P.Value < 1e-4]) ==
            sign(sex_pm_standard$t[sex_pm_standard$P.Value < 1e-4])
)
mean(sex_liv_cell$P.Value[sex_pm_cell$P.Value < 1e-4] < 0.05 &
    sign(sex_liv_cell$t[sex_pm_cell$P.Value < 1e-4]) ==
        sign(sex_pm_cell$t[sex_pm_cell$P.Value < 1e-4]))


smoothScatter(sex_liv$t, sex_pm$t)
smoothScatter(sex_liv_standard$t, sex_pm_standard$t)
smoothScatter(sex_liv_cell$t, sex_pm_cell$t)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
