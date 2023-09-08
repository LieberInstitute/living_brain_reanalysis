library("here")
library("SummarizedExperiment")
library("readxl")
library("jaffelab")
library("edgeR")
library("sessioninfo")

## For reproducing the jitter output
set.seed(20230907)

## Create output directories
dir_plots <- here("plots", "05_check_gene_exprs")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir_rdata <- here("processed-data", "05_check_gene_exprs")
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
cellProps <- read.csv(file.path(dir_rdata, "LBP_burkeDecon.csv"), row.names = 1)
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

## For simplicity later on, I'll use this varaible name
rse_gene$samplingAge <- rse_gene$ageDeath_num.biospecimen

## model
mod <- model.matrix(
    ~ COI + race + sex + diagnosis +
        rin + cellPC +
        mitoRate +
        totalAssignedGene +
        rRNA_rate +
        overallMapRate,
    data = colData(rse_gene)
)
qsvs_cell <- read.csv("qsv_matrices/qSVs_cell.csv", row.names = 1)
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
table(outGene_qsva$adj.P.Val < 0.05)
table(outGene_pc3$adj.P.Val < 0.05)

outGene_order <- outGene[order(outGene$P.Value), ]
outGene_order$mean_dir <- cumsum(sign(outGene_order$t)) / 1:nrow(outGene_order)
plot(outGene_order$mean_dir)

pdf("plots/published_vs_redone_tstat_modelCompare.pdf")
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

pdf("plots/qsva_attentuation.pdf")
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

## model
mod_liv <- model.matrix(
    ~ samplingAge + race + sex + diagnosis +
        rin + recount_qc.star.uniquely_mapped_reads_perc + cellPC +
        recount_qc.aligned_readsperc.chrm +
        recount_qc.aligned_readsperc.chrx +
        recount_qc.bc_auc.all_perc +
        recount_qc.star.perc_of_reads_mapped_to_multiple_loci +
        recount_qc.gene_fc.all_perc +
        recount_qc.junction_avg_coverage,
    data = colData(rse_liv)
)
qsvs_liv_cell <- read.csv("qsv_matrices/qSVs_cell_liv.csv", row.names = 1)
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
    coef = 2,
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_liv$isExprs <- outGene_liv$ensembl_id %in% res$ensembl_id
outGene_liv$adj.P.Val <- NA
outGene_liv$adj.P.Val[outGene_liv$isExprs] <- p.adjust(outGene_liv$P.Value[outGene_liv$isExprs], "fdr")

outGene_qsva_liv <- topTable(
    eBGene_liv_qsva,
    coef = 2,
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_qsva_liv$isExprs <- outGene_qsva_liv$ensembl_id %in% res$ensembl_id
outGene_qsva_liv$adj.P.Val <- NA
outGene_qsva_liv$adj.P.Val[outGene_qsva_liv$isExprs] <- p.adjust(outGene_qsva_liv$P.Value[outGene_qsva_liv$isExprs], "fdr")

table(outGene_liv$adj.P.Val < 0.05)
table(outGene_qsva_liv$adj.P.Val < 0.05)

###########
## age in pm
rse_pm <- rse_gene[, rse_gene$COI == "PM"]

## model
mod_pm <- model.matrix(
    ~ samplingAge + race + sex + diagnosis +
        rin + recount_qc.star.uniquely_mapped_reads_perc + cellPC +
        recount_qc.aligned_readsperc.chrm +
        recount_qc.aligned_readsperc.chrx +
        recount_qc.bc_auc.all_perc +
        recount_qc.star.perc_of_reads_mapped_to_multiple_loci +
        recount_qc.gene_fc.all_perc +
        recount_qc.junction_avg_coverage,
    data = colData(rse_pm)
)
qsvs_pm_cell <- read.csv("qsv_matrices/qSVs_cell_pm.csv", row.names = 1)
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
    coef = 2,
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_pm$isExprs <- outGene_pm$ensembl_id %in% res$ensembl_id
outGene_pm$adj.P.Val <- NA
outGene_pm$adj.P.Val[outGene_pm$isExprs] <- p.adjust(outGene_pm$P.Value[outGene_pm$isExprs], "fdr")

outGene_qsva_pm <- topTable(
    eBGene_pm_qsva,
    coef = 2,
    p.value = 1,
    sort = "none",
    number = nrow(dge)
)
outGene_qsva_pm$isExprs <- outGene_qsva_pm$ensembl_id %in% res$ensembl_id
outGene_qsva_pm$adj.P.Val <- NA
outGene_qsva_pm$adj.P.Val[outGene_qsva_pm$isExprs] <- p.adjust(outGene_qsva_pm$P.Value[outGene_qsva_pm$isExprs], "fdr")

table(outGene_pm$adj.P.Val < 0.05)
table(outGene_qsva_pm$adj.P.Val < 0.05)


plot(
    outGene_pm$logFC[outGene_pm$isExprs],
    outGene_liv$logFC[outGene_pm$isExprs],
    xlab = "Age logFC (PM, no qSVA)",
    ylab = "Age logFC (LIV, no qSVA)",
    pch = 21,
    bg = "grey"
)
cor.test(outGene_pm$logFC[outGene_pm$isExprs], outGene_liv$logFC[outGene_pm$isExprs])
plot(
    outGene_qsva_pm$logFC[outGene_pm$isExprs],
    outGene_qsva_liv$logFC[outGene_pm$isExprs],
    xlab = "Age logFC (PM, qSVA)",
    ylab = "Age logFC (LIV, qSVA)",
    pch = 21,
    bg = "grey"
)
cor.test(outGene_qsva_pm$logFC[outGene_pm$isExprs], outGene_qsva_liv$logFC[outGene_pm$isExprs])

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
