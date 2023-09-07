library("readxl")
library("jaffelab")
library("edgeR")
library("here")
library("sessioninfo")

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
    rpkm(dge, gene.length = rowData(rse_gene)$score, log = TRUE)
pca <- prcomp(t(log_rpkm))
pcaVars <- getPcaVars(pca)

pcHeatmap(rse_gene, pca, plotdir = dir_plots, pthresh = 20)

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

boxplot(rse_gene$mitoRate ~ rse_gene$COI)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
