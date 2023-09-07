library("readxl")
library("jaffelab")
library("RColorBrewer")
library("here")
library("SummarizedExperiment")
library("sessioninfo")

## read in decon
decon_df <- read_excel("decon_model.xlsx", sheet = 1, skip = 3)
colnames(decon_df)[3] <- "CellType"

## make matrix
coefEsts <- as.matrix(decon_df[, 4:13])
rownames(coefEsts) <- ss(decon_df$Gene, "\\.")

### get monorail data for qc
options(recount3_url = "https://neuro-recount-ds.s3.amazonaws.com/recount3")
hp <- available_projects()
rse_gene <- create_rse(hp[hp$project == "LBP", ])
seqlevels(rse_gene, pruning.mode = "coarse") <- paste0("chr", c(1:22, "X", "Y", "M"))
assay(rse_gene, "counts") <- transform_counts(rse_gene)
assay(rse_gene, "raw_counts") <- NULL # drop raw counts
rowData(rse_gene)$ensembl_id <- ss(rownames(rse_gene), "\\.")

exprs <- assays(rse_gene)$counts
exprs <- exprs[rowSums(exprs) > 0, ]
# write.csv(exprs, file = gzfile("tables/gene_counts_LBP.csv.gz")) # for bernie

## add coi info
pheno <- read.csv("phenotype/LBP_merged_phenotype.csv")
rse_gene$COI <- factor(pheno$COI[match(colnames(rse_gene), pheno$specimenID)])

## get RPKM
yExprs <- log2(recount::getRPKM(rse_gene, "bp_length") + 1)
rownames(yExprs) <- rowData(rse_gene)$ensembl_id
yExprs_scaled <- scale(yExprs[rownames(coefEsts), ]) # filter to cell type genes, and scale

## do deconvolution
propEsts <- minfi:::projectCellType(yExprs_scaled, coefEsts)
propEsts_scaled <- prop.table(propEsts, 1)
write.csv(propEsts_scaled, file = "phenotype/LBP_burkeDecon.csv")

##  each variable vs COI
pdf("plots/cellTypes_vs_coi.pdf")
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
