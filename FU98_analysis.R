# FU98 analysis on uppmax

#rRNA rates
setwd("~/Documents/FU98/analysis/FU98_analysis")
rna_qc_uppsala <- read.delim("../uppsala/pipeline_output/RNA_QC/aggregated_metrics.tsv")
rna_qc_stockholm <- read.delim("../stockholm/pipeline_output/RNA_QC/aggregated_metrics.tsv")
rna_qc_uppsala$center <- "Uppsala"
rna_qc_uppsala$pretty.sample <- c("Ua.HeLa-S3.1", "Ua.GM12878.1", "Ua.GM12878.2")
rna_qc_stockholm$center <- "Stockholm"
rna_qc_stockholm$pretty.sample <- c("Sth.GM12878.1", "Sth.HeLa-S3.1", "Sth.GM12878.2", "Sth.HeLa-S3.2")
rna_qc_all <- rbind(rna_qc_uppsala, rna_qc_stockholm)
rna_qc_all$center <- as.factor(rna_qc_all$center)
rna_qc_all.sorted <- rna_qc_all[with(rna_qc_all, order(pretty.sample)), ]

row.names(rna_qc_all.sorted) <- rna_qc_all.sorted$pretty.sample
old <- par()
par(mar=c(10,5,5,5))
barplot (rna_qc_all.sorted$rRNA.rate, col = c("lightgreen", "lightblue")[rna_qc_all.sorted$center], names.arg = row.names(rna_qc_all.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.004), ylab = "rRNA rate")

# Strand specificity
rna_qc_all.sorted$strand.spec <- rna_qc_all.sorted$End.1.Sense / rna_qc_all.sorted$End.2.Sense
rna_qc_all.sorted$strand.spec.log10 <- log10(rna_qc_all.sorted$strand.spec)
barplot(rna_qc_all.sorted$strand.spec.log10, col = c("lightgreen", "lightblue")[rna_qc_all.sorted$center], names.arg = row.names(rna_qc_all.sorted), las = 2, cex.names  = 0.9, ylim = c(-2.5,0.1), ylab = "Strand specificity")
