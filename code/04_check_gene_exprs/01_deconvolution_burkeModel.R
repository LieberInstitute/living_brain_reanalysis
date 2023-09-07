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
dir_plots <- here("plots", "04_check_gene_exprs")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir_rdata <- here("processed-data", "04_check_gene_exprs")
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
