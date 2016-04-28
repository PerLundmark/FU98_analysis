# FU98 analysis
#vars
analysis.directory <- "~/Documents/FU98/analysis/FU98_analysis"
flagstat.directory <- "/Users/per/Documents/FU98/analysis/duplicates/"
longitudinal.flagstat.directory <- "/Users/per/Documents/FU98/analysis/duplicates/"

#Functions
flagstat.mapped.reads <- function(flagstat.file){
  flagstat.df <- read.delim(flagstat.file, sep="", header=F, colClasses = c("character"))
  as.numeric(flagstat.df[5,1])
}

flagstat.mapped.reads.old <- function(flagstat.file){
  flagstat.df <- read.delim(flagstat.file, sep="", header=F, colClasses = c("character"))
  as.numeric(flagstat.df[3,1])
}

flagstat.duplicate.reads <- function(flagstat.file){
  flagstat.df <- read.delim(flagstat.file, sep="", header=F, colClasses = c("character"))
  as.numeric(flagstat.df[4,1])
}

flagstat.duplicate.reads.old <- function(flagstat.file){
  flagstat.df <- read.delim(flagstat.file, sep="", header=F, colClasses = c("character"))
  as.numeric(flagstat.df[2,1])
}

add.flagstat.path <- function(flagstat.filename){
  paste(flagstat.directory,flagstat.filename, sep="")
}

###########################
# Uppsala vs Stockholm 2015
###########################

#rRNA rates
setwd(analysis.directory)
rna_qc_uppsala <- read.delim("../uppsala/pipeline_output/RNA_QC/aggregated_metrics.tsv")
rna_qc_stockholm <- read.delim("../stockholm/pipeline_output/RNA_QC/aggregated_metrics.tsv")
rna_qc_uppsala$center <- "Uppsala"
rna_qc_uppsala$pretty.sample <- c("Ua.HeLa-S3.1", "Ua.GM12878.1", "Ua.GM12878.2")
rna_qc_stockholm$center <- "Stockholm"
rna_qc_stockholm$pretty.sample <- c("Sth.GM12878.1", "Sth.HeLa-S3.2", "Sth.GM12878.2", "Sth.HeLa-S3.1")
rna_qc_all <- rbind(rna_qc_uppsala, rna_qc_stockholm)
rna_qc_all$center <- as.factor(rna_qc_all$center)
rna_qc_all.sorted <- rna_qc_all[with(rna_qc_all, order(pretty.sample)), ]
rna_qc_all.sorted$raw.pcr.dup <- c(0.381183,0.394302,0.297067,0.281671,0.247557,0.217364,0.198524) #not used
flagstat.filenames <- c("ds_1-GM12878-1.bam.flagstat", "ds_2-GM12878-1.bam.flagstat", "ds_1-HeLa-S3.bam.flagstat", "ds_2-HeLa-S3.bam.flagstat", "ds_GM12878-1-RNA.bam.flagstat", "dedup_GM12878-1-RNA-1.bam.flagstat", "ds_HeLs-S3.bam.flagstat" )
rna_qc_all.sorted$flagstat.filename <- flagstat.filenames

row.names(rna_qc_all.sorted) <- rna_qc_all.sorted$pretty.sample

old <- par()
ppi <- 300
png("./output/rRNA_Ua_vs_Sth.png", width=6*ppi, height=6*ppi, res=ppi)
par(mar=c(10,5,5,1))
barplot (rna_qc_all.sorted$rRNA.rate, col = c("lightgreen", "lightblue")[rna_qc_all.sorted$center], names.arg = row.names(rna_qc_all.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.004), ylab = "rRNA rate")
dev.off()

# Strand specificity
rna_qc_all.sorted$strand.spec <- rna_qc_all.sorted$End.1.Sense / rna_qc_all.sorted$End.2.Sense
rna_qc_all.sorted$strand.spec.log10 <- log10(rna_qc_all.sorted$strand.spec)

ppi <- 300
png("./output/Strand_spec_Ua_vs_Sth.png", width=6*ppi, height=6*ppi, res=ppi)
par(mar=c(10,5,5,1))
barplot(rna_qc_all.sorted$strand.spec.log10, col = c("lightgreen", "lightblue")[rna_qc_all.sorted$center], names.arg = row.names(rna_qc_all.sorted), las = 2, cex.names  = 0.9, ylim = c(-2.5,0.1), ylab = "Strand specificity")
dev.off()

# PCR duplicate levels
rna_qc_all.sorted$ds.pcr.dup <- sapply(add.flagstat.path(rna_qc_all.sorted$flagstat.filename), flagstat.duplicate.reads) / sapply(add.flagstat.path(rna_qc_all.sorted$flagstat.filename), flagstat.mapped.reads)
ppi <- 300
png("./output/PCR_duprate_Ua_vs_Sth.png", width=6*ppi, height=6*ppi, res=ppi)
par(mar=c(10,5,5,1))
barplot(rna_qc_all.sorted$ds.pcr.dup, col = c("lightgreen", "lightblue")[rna_qc_all.sorted$center], names.arg = row.names(rna_qc_all.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.4), ylab = "Duplication rate")
dev.off()

# FPKM correlations
library(cummeRbund)

#GM12878
separate.cuff.GM12878 <- readCufflinks("../diff_out_GM12878_masked_separate")
#plot scatterplots of all GM samples together

ppi <- 300
png("./output/Scatter_GM12878.png", width=8*ppi, height=8*ppi, res=ppi)
csScatterMatrix(genes(separate.cuff.GM12878))
dev.off()
#sc.plot1 <- csScatter(genes(separate.cuff), 'Uppsala_GM12878_1', 'Stockholm_GM12878_1') + annotate("text", x = 4, y = 1e6, label = "Pearson correlation: X", size = 10) + theme(axis.title = element_text(size=20))
#sc.plot1 + annotate("text", x = 4, y = 1e6, label = "Spearman correlation: X", size = 10)

#Table of Spearman correlations
separate.fpkm.GM12878 <- fpkmMatrix(genes(separate.cuff.GM12878))
spearman.corr.GM12878 <- cor(separate.fpkm.GM12878, method = "spearman")
write.table(spearman.corr.GM12878, file="./output/spearman.corr.GM12878.tsv", sep="\t")

#Pearson correlations
cor(log2(separate.fpkm.GM12878 + 1), method = "pearson")
cor(log2(separate.fpkm.GM12878 + .01), method = "pearson")
cor(log2(separate.fpkm.GM12878 + .001), method = "pearson")
cor(log10(separate.fpkm.GM12878 + 1), method = "pearson")
cor(log10(separate.fpkm.GM12878 + .01), method = "pearson")
cor(log10(separate.fpkm.GM12878 + .001), method = "pearson")



#HeLa-S3
separate.cuff.HeLa <- readCufflinks("../diff_out_HeLa_masked_separate")
#plot scatterplots of all GM samples together
ppi <- 300
png("./output/Scatter_HeLa.png", width=8*ppi, height=8*ppi, res=ppi)
csScatterMatrix(genes(separate.cuff.HeLa))
dev.off()
#sc.plot1 <- csScatter(genes(separate.cuff), 'Uppsala_GM12878_1', 'Stockholm_GM12878_1') + annotate("text", x = 4, y = 1e6, label = "Pearson correlation: X", size = 10) + theme(axis.title = element_text(size=20))
#sc.plot1 + annotate("text", x = 4, y = 1e6, label = "Spearman correlation: X", size = 10)

#Table of Spearman correlations
separate.fpkm.HeLa <- fpkmMatrix(genes(separate.cuff.HeLa))
spearman.corr.HeLa <- cor(separate.fpkm.HeLa, method = "spearman")
write.table(spearman.corr.HeLa, file="./output/spearman.corr.HeLa.tsv", sep="\t")

#pearson correlations
cor(log2(separate.fpkm.HeLa + 1), method = "pearson")
cor(log2(separate.fpkm.HeLa + .01), method = "pearson")
cor(log2(separate.fpkm.HeLa + .001), method = "pearson")
cor(log10(separate.fpkm.HeLa + 1), method = "pearson")
cor(log10(separate.fpkm.HeLa + .01), method = "pearson")
cor(log10(separate.fpkm.HeLa + .001), method = "pearson")



#########################
### Longitudinal analysis
#########################

##rRNA rates
rna_qc_2013 <- read.delim("../FU45/RNA_QC/aggregated_metrics.tsv")
rna_qc_2014 <- read.delim("../uppsala_2014/pipeline_output/RNA_QC/aggregated_metrics.tsv")
rna_qc_2013$center <- "Uppsala-2013"
rna_qc_2014$center <- "Uppsala-2014"
rna_qc_2013$pretty.sample <- c("Ua-2013.HeLa-S3.1ug.1", "Ua-2013.GM12878.100ng.2", "Ua-2013.GM12878.1ug.1", "Ua-2013.GM12878.100ng.1", "Ua-2013.HeLa-S3.1ug.2", "Ua-2013.GM12878.1ug.2")
rna_qc_2014$pretty.sample <- c("Ua-2014.HeLa-S3.2", "Ua-2014.HeLa-S3.1");

rna_qc_all_old_vs_new <- rbind(rna_qc_uppsala, rna_qc_2013)
rna_qc_all_2015_2014_2013 <- rbind(rna_qc_all_old_vs_new, rna_qc_2014)
#rna_qc_all_old_vs_new$center <- as.factor(rna_qc_all_old_vs_new$center)
#rna_qc_all_old_vs_new.sorted <- rna_qc_all_old_vs_new[with(rna_qc_all_old_vs_new, order(pretty.sample)), ]
rna_qc_all_2015_2014_2013$center <- as.factor(rna_qc_all_2015_2014_2013$center)
rna_qc_all_2015_2014_2013.sorted <- rna_qc_all_2015_2014_2013[with(rna_qc_all_2015_2014_2013, order(pretty.sample)), ]
flagstat.filenames.longitudinal <- c(
"/Users/per/Documents/FU98/analysis/FU45/flagstat/GM12878_Str-polyA_100ng-1.flagstat.txt",
"/Users/per/Documents/FU98/analysis/FU45/flagstat/GM12878_Str-polyA_100ng-2.flagstat.txt",
"/Users/per/Documents/FU98/analysis/FU45/flagstat/GM12878_Str-polyA_1ug-1.flagstat.txt", 
"/Users/per/Documents/FU98/analysis/FU45/flagstat/GM12878_Str-polyA_1ug-2.flagstat.txt",
"/Users/per/Documents/FU98/analysis/FU45/flagstat/HeLa_Str-polyA_1ug-1.flagstat.txt",
"/Users/per/Documents/FU98/analysis/FU45/flagstat/HeLa_Str-polyA_1ug-2.flagstat.txt",
"/Users/per/Documents/FU98/analysis/duplicates_longitudinal/ds_markdup_2014_HeLa-1.bam.flagstat",
"/Users/per/Documents/FU98/analysis/duplicates_longitudinal/ds_markdup_2014_HeLa-2.bam.flagstat",
"/Users/per/Documents/FU98/analysis/duplicates_longitudinal/ds_markdup_2015_GM12878-1.bam.flagstat",
"/Users/per/Documents/FU98/analysis/duplicates_longitudinal/ds_markdup_2015_GM12878-2.bam.flagstat",
"/Users/per/Documents/FU98/analysis/duplicates_longitudinal/ds_markdup_2015_HeLa.bam.flagstat"
)
rna_qc_all_2015_2014_2013.sorted$flagstat.filename <- flagstat.filenames.longitudinal

#row.names(rna_qc_all_old_vs_new.sorted) <- rna_qc_all_old_vs_new.sorted$pretty.sample
row.names(rna_qc_all_2015_2014_2013.sorted) <-  rna_qc_all_2015_2014_2013.sorted$pretty.sample

old <- par()

ppi <- 300
png("./output/rRNA_longitudinal.png", width=6*ppi, height=6*ppi, res=ppi)
par(mar=c(12,5,5,1))
#barplot (rna_qc_all_old_vs_new.sorted$rRNA.rate, col = c("lightgreen", "lightblue")[rna_qc_all_old_vs_new.sorted$center], names.arg = row.names(rna_qc_all_old_vs_new.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.004), ylab = "rRNA rate")
barplot (rna_qc_all_2015_2014_2013.sorted$rRNA.rate, col = c("lightgreen", "lightblue", "white")[rna_qc_all_2015_2014_2013.sorted$center], names.arg = row.names(rna_qc_all_2015_2014_2013.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.004), ylab = "rRNA rate")
dev.off()

# Strand specificity
#rna_qc_all_old_vs_new.sorted$strand.spec <- rna_qc_all_old_vs_new.sorted$End.1.Sense / rna_qc_all_old_vs_new.sorted$End.2.Sense
#rna_qc_all_old_vs_new.sorted$strand.spec.log10 <- log10(rna_qc_all_old_vs_new.sorted$strand.spec)
#barplot(rna_qc_all_old_vs_new.sorted$strand.spec.log10, col = c("lightgreen", "lightblue")[rna_qc_all_old_vs_new.sorted$center], names.arg = row.names(rna_qc_all_old_vs_new.sorted), las = 2, cex.names  = 0.9, ylim = c(-2.5,0.1), ylab = "Strand specificity")

rna_qc_all_2015_2014_2013.sorted$strand.spec <- rna_qc_all_2015_2014_2013.sorted$End.1.Sense / rna_qc_all_2015_2014_2013.sorted$End.2.Sense
rna_qc_all_2015_2014_2013.sorted$strand.spec.log10 <- log10(rna_qc_all_2015_2014_2013.sorted$strand.spec)

ppi <- 300
png("./output/Strand_spec_longitudinal.png", width=6*ppi, height=6*ppi, res=ppi)
par(mar=c(12,5,5,1))
barplot(rna_qc_all_2015_2014_2013.sorted$strand.spec.log10, col = c("lightgreen", "lightblue", "white")[rna_qc_all_2015_2014_2013.sorted$center], names.arg = row.names(rna_qc_all_2015_2014_2013.sorted), las = 2, cex.names  = 0.9, ylim = c(-2.5,0.1), ylab = "Strand specificity")
dev.off()

#PCR duplicate levels   work in progress!####
rna_qc_all_2015_2014_2013.sorted$ds.pcr.dup[1:6] <- sapply(rna_qc_all_2015_2014_2013.sorted$flagstat.filename[1:6], flagstat.duplicate.reads.old) / sapply(rna_qc_all_2015_2014_2013.sorted$flagstat.filename[1:6], flagstat.mapped.reads.old)
rna_qc_all_2015_2014_2013.sorted$ds.pcr.dup[7:11] <- sapply(rna_qc_all_2015_2014_2013.sorted$flagstat.filename[7:11], flagstat.duplicate.reads) / sapply(rna_qc_all_2015_2014_2013.sorted$flagstat.filename[7:11], flagstat.mapped.reads)

ppi <- 300
png("./output/PCR_duprate_longitudinal.png", width=6*ppi, height=6*ppi, res=ppi)
par(mar=c(10,5,5,1))
barplot(rna_qc_all_2015_2014_2013.sorted$ds.pcr.dup, col = c("lightgreen", "lightblue")[rna_qc_all_2015_2014_2013.sorted$center], names.arg = row.names(rna_qc_all_2015_2014_2013.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.4), ylab = "Duplication rate")
dev.off()
####

## Correlations / scatterplots
old_new.cuff.GM12878 <- readCufflinks("../diff_out_GM12878_masked_separate_old_vs_new")
old_new.cuff.HeLa <- readCufflinks("../diff_out_HeLa_masked_separate_old_vs_new")
cuff.2015.2014.2013.HeLa <- readCufflinks("../diff_out_HeLa_masked_2013_2014_2015/")


#GM12878
ppi <- 300
png("./output/Scatter_GM12878_longitudinal.png", width=10*ppi, height=10*ppi, res=ppi)
csScatterMatrix(genes(old_new.cuff.GM12878))
dev.off()

separate.fpkm.GM12878.old.new <- fpkmMatrix(genes(old_new.cuff.GM12878))
spearman.corr.GM12878.longitudinal <- cor(separate.fpkm.GM12878.old.new, method = "spearman")
write.table(spearman.corr.GM12878.longitudinal, file="./output/spearman.corr.GM12878.longitudinal.tsv", sep="\t")


#pearson correlations
cor(log2(separate.fpkm.GM12878.old.new + 1), method = "pearson")
cor(log2(separate.fpkm.GM12878.old.new + .01), method = "pearson")
cor(log2(separate.fpkm.GM12878.old.new + .001), method = "pearson")
cor(log10(separate.fpkm.GM12878.old.new + 1), method = "pearson")
cor(log10(separate.fpkm.GM12878.old.new + .01), method = "pearson")
cor(log10(separate.fpkm.GM12878.old.new + .001), method = "pearson")



#HeLa-S3
#csScatterMatrix(genes(old_new.cuff.HeLa))
#separate.fpkm.HeLa.old.new <- fpkmMatrix(genes(old_new.cuff.HeLa))
#cor(separate.fpkm.HeLa.old.new, method = "spearman")

ppi <- 300
png("./output/Scatter_HeLa_longitudinal.png", width=10*ppi, height=10*ppi, res=ppi)
csScatterMatrix(genes(cuff.2015.2014.2013.HeLa))
dev.off()
separate.fpkm.HeLa.2015.2014.2013 <- fpkmMatrix(genes(cuff.2015.2014.2013.HeLa))
spearman.corr.HeLa.longitudinal <- cor(separate.fpkm.HeLa.2015.2014.2013, method = "spearman")
write.table(spearman.corr.HeLa.longitudinal, file="./output/spearman.corr.HeLa.longitudinal.tsv", sep="\t")

#pearson correlations
cor(log2(separate.fpkm.HeLa.2015.2014.2013 + 1), method = "pearson")
cor(log2(separate.fpkm.HeLa.2015.2014.2013 + .01), method = "pearson")
cor(log2(separate.fpkm.HeLa.2015.2014.2013 + .001), method = "pearson")
cor(log10(separate.fpkm.HeLa.2015.2014.2013 + 1), method = "pearson")
cor(log10(separate.fpkm.HeLa.2015.2014.2013 + .01), method = "pearson")
cor(log10(separate.fpkm.HeLa.2015.2014.2013 + .001), method = "pearson")
