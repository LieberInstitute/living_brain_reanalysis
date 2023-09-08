library("qsvaR")
library("SummarizedExperiment")
library("parallel")
library("jaffelab")
library("limma")
library("rtracklayer")
library("recount3")
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

dir.create("qsv_matrices")
dir.create("phenotype")

data_path = "/shared-data/research/genomics/datasets/synapse/LBP/"

## read in transcripts
tpm = read.delim(file.path(data_path, "recount_salmon/all_salmon_tpms.tsv"),
	row.names=1)
tpm$Length = NULL

## annotation info
gtf = import("/shared-data/research/genomics/datasets/recount3/ref/hg38/gtf/genes.transcripts.gtf")
anno = mcols(gtf)[,c("transcript_id", "transcript_type", "gene_name", "gene_id")]
anno = anno[match(rownames(tpm), anno$transcript_id),]

# ## get pheno for postmortem vs living
# biospec = read.csv(file.path(data_path, "LBP_biospecimen_metadata.csv"))
# assay = read.csv(file.path(data_path, "LBP_assay_RNAseq_metadata.csv"))
# indiv = read.csv(file.path(data_path, "LBP_individual_metadata.csv"))
# pheno = dplyr::left_join(assay, biospec)
# pheno = dplyr::left_join(pheno, indiv)
# write.csv(pheno, file = "phenotype/LBP_merged_phenotype.csv",row.names=FALSE)
pheno = read.csv("phenotype/LBP_merged_phenotype.csv")
pheno$COI = factor(ifelse(pheno$isPostMortem, "PM", "LIV"))

## put tpm in same order
tpm = tpm[,pheno$specimenID]

## read in degradation tpm
tpm_degrade = read.delim("/shared-data/research/genomics/datasets/libd/SRP108559/SRP108559.TPM.pasted.tsv",
	row.names=1)
tpm_degrade = tpm_degrade[,-1] # drop length
hp = available_projects()
rse_degrade = create_rse(hp[hp$project == 'SRP108559',])
rse_degrade = expand_sra_attributes(rse_degrade)
rse_degrade = rse_degrade[,rse_degrade$sra_attribute.library_type == "RiboZero"]
pheno_degrade = colData(rse_degrade)
pheno_degrade$BrNum = factor(pheno_degrade$sra_attribute.brain_number)
pheno_degrade$degrade_time = as.numeric(pheno_degrade$sra_attribute.degradation_time)

pheno_degrade = pheno_degrade[order(pheno_degrade$BrNum, pheno_degrade$degrade_time),]
tpm_degrade = tpm_degrade[,rownames(pheno_degrade)]

### recalculate degradation stas for all isoforms, not just 45k in objects
mod_degrade = model.matrix(~degrade_time + BrNum, data = pheno_degrade)
fit_degrade = lmFit(log2(tpm_degrade+1), mod_degrade)
degradation_tstats = topTable(eBayes(fit_degrade), coef=2,
	sort = "none", n = nrow(tpm_degrade), genelist = anno)
degradation_tstats = degradation_tstats[order(degradation_tstats$P.Value),]


## plots for top degradation t-stats
# filter to rows
cell_rows = rownames(tpm) %in% transcripts$cell_component
standard_rows = rownames(tpm) %in% transcripts$standard

top10_rows = rownames(tpm) %in% rownames(degradation_tstats)[1:10]
top25_rows = rownames(tpm) %in% rownames(degradation_tstats)[1:25]
top50_rows = rownames(tpm) %in% rownames(degradation_tstats)[1:50]
top100_rows = rownames(tpm) %in% rownames(degradation_tstats)[1:100]
top200_rows = rownames(tpm) %in% rownames(degradation_tstats)[1:200]

## make list of TPMs
tpm_list = list(tpm_standard = tpm[standard_rows,],
				tpm_standard_pm = tpm[standard_rows,pheno$COI == "PM"],
				tpm_standard_liv = tpm[standard_rows,pheno$COI == "LIV"],
				tpm_cell = tpm[cell_rows,],
				tpm_cell_pm = tpm[cell_rows,pheno$COI == "PM"],
				tpm_cell_liv = tpm[cell_rows,pheno$COI == "LIV"],
				tpm_top10 = tpm[top10_rows,],
				tpm_top25 = tpm[top25_rows,],
				tpm_top50 = tpm[top50_rows,],
				tpm_top100 = tpm[top100_rows,],
				tpm_top200 = tpm[top200_rows,])

set.seed(743)
qsva_list = mclapply(tpm_list, function(x) {
	expr = log2(x + 1)
	m = matrix(1, nrow = ncol(expr), ncol = 1)
	k = sva::num.sv(expr, m)
	pca = prcomp(t(expr))
	qsvs = pca$x[,1:k]
	var_expl = getPcaVars(pca)[1:k]
	list(qSVs = qsvs, varExpl = var_expl)

},mc.cores=6)

save(qsva_list, file = "qsv_matrices/qsva_objects.RData")


#############
## load objects
load("qsv_matrices/qsva_objects.RData")

## write out as text
for(i in seq(along=qsva_list)) {
	write.csv(qsva_list[[i]]$qSVs, file = paste0("qsv_matrices/qSVs_",
			gsub("tpm_", "", names(qsva_list)[i]), ".csv"))
}

# associations to PM vs LIV
coi_assoc = lapply(qsva_list[c(1,4,7:11)], function(y) {
	round(-log10(apply(t(t(y$qSVs)), 2,
		function(x) t.test(x ~ pheno$COI)$p.value)),1)
})
coi_assoc

anno[match(rownames(degradation_tstats)[1:10], anno$transcript_id),]

boxplot(qsva_list$tpm_top10$qSVs ~ pheno$COI,
	ylab = paste0("qSV1 (top10 Tx: ", qsva_list$tpm_top10$varExpl[1],
		"% Var Expl)"),xlab="")


### check qSV vs COI
apply(qsva_list$tpm_standard$qSVs, 2, function(x) t.test(x ~ pheno$COI)$p.value)
apply(qsva_list$tpm_cell$qSVs, 2, function(x) t.test(x ~ pheno$COI)$p.value)

## top qSV is very correlated
cor(qsva_list$tpm_standard$qSVs[,1], qsva_list$tpm_cell$qSVs[,1])

qsva_list$tpm_standard$varExpl


boxplot(qsva_list$tpm_standard$qSVs[,1] ~ pheno$COI,
	ylab = paste0("qSV1 (standard: ", qsva_list$tpm_standard$varExpl[1],
		"% Var Expl)"),xlab="")
boxplot(qsva_list$tpm_cell[,1] ~ pheno$COI,
	ylab = "qSV1 (cell)",xlab="")


### read in rnaseq metrics
metrics = read.csv("phenotype/LBP_merged_seqmetrics.csv")
metrics = metrics[match(pheno$specimenID, metrics$external_id),]
plot(qsva_list$tpm_top10$qSVs, metrics$recount_qc.intron_sum_perc)
plot(qsva_list$tpm_top10$qSVs, metrics$recount_qc.exon_fc.unique_perc)
plot(qsva_list$tpm_top10$qSVs, metrics$"recount_qc.bc_auc.all_reads_annotated_bases")

###########
## de checks
pheno$race[pheno$race == "White "] = "White"
pheno$race[pheno$race == ""] = "Unknown"
pheno$diagnosis[pheno$diagnosis == "Parkinson's disease "] = "Parkinson's disease"
pheno$samplingAge[pheno$samplingAge == "89+"] = 89
pheno$samplingAge = as.numeric(pheno$samplingAge)
rownames(pheno) = pheno$specimenID
pheno$diagnosis = factor(pheno$diagnosis)
pheno$diagnosis = relevel(pheno$diagnosis, "control")

pheno = cbind(pheno, metrics)

## model
mod = model.matrix(~COI + race + sex + diagnosis +
		rin + cellPC + recount_qc.star.uniquely_mapped_reads_perc +
		recount_qc.aligned_readsperc.chrm +
		recount_qc.aligned_readsperc.chrx +
		recount_qc.bc_auc.all_perc +
		recount_qc.star.perc_of_reads_mapped_to_multiple_loci +
		recount_qc.gene_fc.all_perc +
		recount_qc.junction_avg_coverage,
		data=pheno)

## modeling
f = lmFit(log2(tpm+1), mod)
out =  topTable(eBayes(f), coef=2, sort="none", n = nrow(tpm), genelist = anno)

## against degradation
out_top10000 = out[order(out$P.Value)[1:10000],]
degradation_tstats_top10000 = degradation_tstats[
	match(ss(out_top10000$transcript_id, "\\."),
		ss(rownames(degradation_tstats), "\\.")),]
plot(out_top10000$t, degradation_tstats_top10000$t)
tt_10000 = table(`PM>LIV` = out_top10000$t > 0, `Degrade>0` = degradation_tstats_top10000$t > 0)
tt_10000
getOR(tt_10000) # 13.6


#####
f_standard = lmFit(log2(tpm+1), cbind(mod,qsva_list$tpm_standard$qSVs))
out_standard =  topTable(eBayes(f_standard), coef=2, sort="none",
		n = nrow(tpm), genelist = anno)

f_cell = lmFit(log2(tpm+1), cbind(mod,qsva_list$tpm_cell$qSVs))
out_cell =  topTable(eBayes(f_cell), coef=2,
		sort="none", n = nrow(tpm), genelist = anno)

cor(out_cell$t, out$t)
smoothScatter(out_cell$t, out$t)
out_cell[order(out_cell$P.Value)[1:25],c("gene_name", "t", "P.Value")]

f_top100 = lmFit(log2(tpm+1), cbind(mod,qsva_list$tpm_top100$qSVs))
out_top100 =  topTable(eBayes(f_top100), coef=2, sort="none", n = nrow(tpm),
	genelist = anno)

table(out$adj.P.Val < 0.05)
table(out_standard$adj.P.Val < 0.05)
table(out_cell$adj.P.Val < 0.05)
table(out_top100$adj.P.Val < 0.05)

degradation_tstats_out = degradation_tstats[
		match(out$transcript_id, rownames(degradation_tstats)),]

smoothScatter(out$t, degradation_tstats_out$t)

sumt = degradation_tstats_out$t + out$t
# cir_ind = c(which.min(sumt), which.max(sumt))

pdf("plots/Dequal_txLevel.pdf")
par(mar=c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main = 2)
plot(out$t, degradation_tstats_out$t,pch=21,bg="grey",
	xlab = "PM vs LIV (t-stat)", ylab = "Degradation (t-stat)")
# points(out$t[cir_ind], degradation_tstats_out$t[cir_ind],
# 	pch = 1, col = c("red","blue"), cex=1.5)
legend("topleft", paste0("Cor=",
	signif(cor(out$t, degradation_tstats_out$t,use="comp"),3)),
	cex=1.5)
abline(h=0,v=0,lty=2)
dev.off()


# ######################
#### example plots #####
plot_ind = c(order(sumt)[1:25], order(sumt, decreasing=TRUE)[1:25])
degradation_tstats_sumt = degradation_tstats_out[plot_ind,]
tpm_degrade_order = as.matrix(tpm_degrade[rownames(degradation_tstats_sumt),])
bIndexes = splitit(pheno_degrade$BrNum)

pdf("plots/degradation_transcripts_top.pdf")
par(mar=c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main = 2)
palette(brewer.pal(5,"Set1"))
for(i in 1:nrow(tpm_degrade_order)) {
	plot(log2(tpm_degrade_order[i,]+1) ~ pheno_degrade$degrade_time,
		pch = 21, bg = pheno_degrade$BrNum,cex=2,
		xlab = "Degradation Time (min)", ylab = "log2(TPM+1)",
		main = paste0(rownames(degradation_tstats_sumt)[i], "\n",
					degradation_tstats_sumt$gene_name[i], " (",
					degradation_tstats_sumt$transcript_type[i], ")"))
	for(j in seq(along=bIndexes)) {
		lines(log2(tpm_degrade_order[i,]+1) ~ pheno_degrade$degrade_time,
				lwd=2,col=j, subset = bIndexes[[j]])
	}
	lc = ifelse(degradation_tstats_sumt$t[i] > 0, "topleft", "topright")
	legend(lc, paste0("p=", signif(degradation_tstats_sumt$P.Value[i],3)), cex=1.7)
}
dev.off()

tpm_lbp_order = as.matrix(tpm[rownames(tpm_degrade_order),])
pdf("plots/degradation_transcripts_lbp.pdf")
par(mar=c(5,6,4,2), cex.axis=2,cex.lab=2,cex.main = 2)
palette(brewer.pal(5,"Dark2"))
for(i in 1:nrow(tpm_degrade_order)) {
	boxplot(log2(tpm_lbp_order[i,]+1) ~ pheno$COI,
		xlab = "", ylab = "log2(TPM+1)",
		ylim = range(log2(tpm_lbp_order[i,]+1)), outline = FALSE,
		main = paste0(rownames(degradation_tstats_sumt)[i], "\n",
					degradation_tstats_sumt$gene_name[i], " (",
					degradation_tstats_sumt$transcript_type[i], ")"))

	points(log2(tpm_lbp_order[i,]+1) ~ jitter(as.numeric(pheno$COI),amount=0.15),
		pch = 21, bg = pheno$COI)

}
dev.off()

### prediction modeling
X = as.data.frame(t(tpm_degrade[rownames(degradation_tstats)[1:10],]))
X_new = as.data.frame(t(tpm[colnames(X),]))
X$Time = pheno_degrade$degrade_time

f_pred = lm(Time ~ . , data = X)

pred_time = predict(f_pred, X_new)

pdf("plots/degradation_prediction_lbp.pdf")
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2,cex.main = 2)
palette(brewer.pal(5,"Dark2"))
boxplot(pred_time ~ pheno$COI,
	xlab = "", ylab = "Predicted Degradation Time (Min)",
	ylim = range(pred_time), outline = FALSE)
points(pred_time ~ jitter(as.numeric(pheno$COI),amount=0.15),
	pch = 21, bg = pheno$COI)
legend("topleft", paste0("p=",
	signif(summary(lm(pred_time ~ pheno$COI))$coef[2,4],3)),cex=2)
dev.off()

##########
## age + sex + dx?
##########

## in living
pheno_liv = pheno[pheno$COI == "LIV",]
pheno_liv$diagnosis = droplevels(pheno_liv$diagnosis )
mod_liv = model.matrix(~samplingAge + sex+ diagnosis + race +   rnaBatch + libraryBatch,
					data = pheno_liv)
tpm_liv = tpm[,rownames(mod_liv)]
standard_qsv_liv = qsva_list$tpm_standard_liv$qSVs[rownames(mod_liv),]
cell_qsv_liv = qsva_list$tpm_cell_liv$qSVs[rownames(mod_liv),]

f_liv = lmFit(log2(tpm_liv+1), mod_liv)
f_liv_standard = lmFit(log2(tpm_liv+1), cbind(mod_liv,standard_qsv_liv))
f_liv_cell = lmFit(log2(tpm_liv+1), cbind(mod_liv,cell_qsv_liv))

## in pm
pheno_pm = pheno[pheno$COI == "PM",]
pheno_pm$diagnosis = droplevels(pheno_pm$diagnosis )
mod_pm = model.matrix(~samplingAge + sex + diagnosis + race + rnaBatch + libraryBatch,
					data = pheno_pm)
tpm_pm = tpm[,rownames(mod_pm)]
standard_qsv_pm = qsva_list$tpm_standard_pm$qSVs[rownames(mod_pm),]
cell_qsv_pm = qsva_list$tpm_cell_pm$qSVs[rownames(mod_pm),]

f_pm = lmFit(log2(tpm_pm+1), mod_pm)
f_pm_standard = lmFit(log2(tpm_pm+1), cbind(mod_pm,standard_qsv_pm))
f_pm_cell = lmFit(log2(tpm_pm+1), cbind(mod_pm,cell_qsv_pm))

#############
##### age ###
age_liv = topTable(eBayes(f_liv), coef=2, sort="none", n = nrow(tpm))
age_liv_standard = topTable(eBayes(f_liv_standard),
	coef=2, sort="none", n = nrow(tpm))
age_liv_cell = topTable(eBayes(f_liv_cell),
	coef=2, sort="none", n = nrow(tpm))
age_pm = topTable(eBayes(f_pm), coef=2, sort="none", n = nrow(tpm))
age_pm_standard = topTable(eBayes(f_pm_standard),
	coef=2, sort="none", n = nrow(tpm))
age_pm_cell = topTable(eBayes(f_pm_cell),
	coef=2, sort="none", n = nrow(tpm))

## liv discovery
mean(age_pm$P.Value[age_liv$P.Value < 1e-4] < 0.05 &
		sign(age_pm$t[age_liv$P.Value < 1e-4]) ==
			sign(age_liv$t[age_liv$P.Value < 1e-4]))
mean(age_pm_standard$P.Value[age_liv_standard$P.Value < 1e-4] < 0.05 &
		sign(age_pm_standard$t[age_liv_standard$P.Value < 1e-4]) ==
			sign(age_liv_standard$t[age_liv_standard$P.Value < 1e-4]))
mean(age_pm_cell$P.Value[age_liv_cell$P.Value < 1e-4] < 0.05 &
		sign(age_pm_cell$t[age_liv_cell$P.Value < 1e-4]) ==
			sign(age_liv_cell$t[age_liv_cell$P.Value < 1e-4]))

## pm discovery
mean(age_liv$P.Value[age_pm$P.Value < 1e-4] < 0.05 &
		sign(age_liv$t[age_pm$P.Value < 1e-4]) ==
			sign(age_pm$t[age_pm$P.Value < 1e-4]))
mean(age_liv_standard$P.Value[age_pm_standard$P.Value < 1e-4] < 0.05 &
		sign(age_liv_standard$t[age_pm_standard$P.Value < 1e-4]) ==
			sign(age_pm_standard$t[age_pm_standard$P.Value < 1e-4]))
mean(age_liv_cell$P.Value[age_pm_cell$P.Value < 1e-4] < 0.05 &
		sign(age_liv_cell$t[age_pm_cell$P.Value < 1e-4]) ==
			sign(age_pm_cell$t[age_pm_cell$P.Value < 1e-4]))


smoothScatter(age_liv$t, age_pm$t)
smoothScatter(age_liv_standard$t, age_pm_standard$t)
smoothScatter(age_liv_cell$t, age_pm_cell$t)

#############
##### sex ###
sex_liv = topTable(eBayes(f_liv), coef=9, sort="none", n = nrow(tpm))
sex_liv_standard = topTable(eBayes(f_liv_standard),
	coef=9, sort="none", n = nrow(tpm))
sex_liv_cell = topTable(eBayes(f_liv_cell),
	coef=9, sort="none", n = nrow(tpm))
sex_pm = topTable(eBayes(f_pm), coef=9, sort="none", n = nrow(tpm))
sex_pm_standard = topTable(eBayes(f_pm_standard),
	coef=9, sort="none", n = nrow(tpm))
sex_pm_cell = topTable(eBayes(f_pm_cell),
	coef=9, sort="none", n = nrow(tpm))

## liv discovery
mean(sex_pm$P.Value[sex_liv$P.Value < 1e-4] < 0.05 &
		sign(sex_pm$t[sex_liv$P.Value < 1e-4]) ==
			sign(sex_liv$t[sex_liv$P.Value < 1e-4]))
mean(sex_pm_standard$P.Value[sex_liv_standard$P.Value < 1e-4] < 0.05 &
		sign(sex_pm_standard$t[sex_liv_standard$P.Value < 1e-4]) ==
			sign(sex_liv_standard$t[sex_liv_standard$P.Value < 1e-4]))
mean(sex_pm_cell$P.Value[sex_liv_cell$P.Value < 1e-4] < 0.05 &
		sign(sex_pm_cell$t[sex_liv_cell$P.Value < 1e-4]) ==
			sign(sex_liv_cell$t[sex_liv_cell$P.Value < 1e-4]))

## pm discovery
mean(sex_liv$P.Value[sex_pm$P.Value < 1e-4] < 0.05 &
		sign(sex_liv$t[sex_pm$P.Value < 1e-4]) ==
			sign(sex_pm$t[sex_pm$P.Value < 1e-4]))
mean(sex_liv_standard$P.Value[sex_pm_standard$P.Value < 1e-4] < 0.05 &
		sign(sex_liv_standard$t[sex_pm_standard$P.Value < 1e-4]) ==
			sign(sex_pm_standard$t[sex_pm_standard$P.Value < 1e-4]))
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
