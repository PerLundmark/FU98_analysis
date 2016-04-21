# FU98 analysis
#vars
analysis.directory <- "~/Documents/FU98/analysis/FU98_analysis"
flagstat.directory <- "/Users/per/Documents/FU98/analysis/duplicates/"

#Functions
flagstat.mapped.reads <- function(flagstat.file){
  flagstat.df <- read.delim(flagstat.file, sep="", header=F, colClasses = c("character"))
  as.numeric(flagstat.df[5,1])
}

flagstat.duplicate.reads <- function(flagstat.file){
  flagstat.df <- read.delim(flagstat.file, sep="", header=F, colClasses = c("character"))
  as.numeric(flagstat.df[4,1])
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
rna_qc_uppsala$center <- "../uppsala/pipeline_output/RNA_QC/aggregated_metrics.tsv"
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
par(mar=c(10,5,5,5))
barplot (rna_qc_all.sorted$rRNA.rate, col = c("lightgreen", "lightblue")[rna_qc_all.sorted$center], names.arg = row.names(rna_qc_all.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.004), ylab = "rRNA rate")

# Strand specificity
rna_qc_all.sorted$strand.spec <- rna_qc_all.sorted$End.1.Sense / rna_qc_all.sorted$End.2.Sense
rna_qc_all.sorted$strand.spec.log10 <- log10(rna_qc_all.sorted$strand.spec)
barplot(rna_qc_all.sorted$strand.spec.log10, col = c("lightgreen", "lightblue")[rna_qc_all.sorted$center], names.arg = row.names(rna_qc_all.sorted), las = 2, cex.names  = 0.9, ylim = c(-2.5,0.1), ylab = "Strand specificity")

# PCR duplicate levels
rna_qc_all.sorted$ds.pcr.dup <- sapply(add.flagstat.path(rna_qc_all.sorted$flagstat.filename), flagstat.duplicate.reads) / sapply(add.flagstat.path(rna_qc_all.sorted$flagstat.filename), flagstat.mapped.reads)

# FPKM correlations
library(cummeRbund)

#GM12878
separate.cuff.GM12878 <- readCufflinks("../diff_out_GM12878_masked_separate")
#plot scatterplots of all GM samples together
csScatterMatrix(genes(separate.cuff.GM12878))
#sc.plot1 <- csScatter(genes(separate.cuff), 'Uppsala_GM12878_1', 'Stockholm_GM12878_1') + annotate("text", x = 4, y = 1e6, label = "Pearson correlation: X", size = 10) + theme(axis.title = element_text(size=20))
#sc.plot1 + annotate("text", x = 4, y = 1e6, label = "Spearman correlation: X", size = 10)

#Table of Spearman correlations
separate.fpkm.GM12878 <- fpkmMatrix(genes(separate.cuff.GM12878))
cor(separate.fpkm.GM12878, method = "spearman")

#HeLa-S3
separate.cuff.HeLa <- readCufflinks("../diff_out_HeLa_masked_separate")
#plot scatterplots of all GM samples together
csScatterMatrix(genes(separate.cuff.HeLa))
#sc.plot1 <- csScatter(genes(separate.cuff), 'Uppsala_GM12878_1', 'Stockholm_GM12878_1') + annotate("text", x = 4, y = 1e6, label = "Pearson correlation: X", size = 10) + theme(axis.title = element_text(size=20))
#sc.plot1 + annotate("text", x = 4, y = 1e6, label = "Spearman correlation: X", size = 10)

#Table of Spearman correlations
separate.fpkm.HeLa <- fpkmMatrix(genes(separate.cuff.HeLa))
cor(separate.fpkm.HeLa, method = "spearman")

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

#from primary above, saving in case its useful later
#rna_qc_all.sorted$raw.pcr.dup <- c(0.381183,0.394302,0.297067,0.281671,0.247557,0.217364,0.198524)
#flagstat.filenames <- c("ds_1-GM12878-1.bam.flagstat", "ds_2-GM12878-1.bam.flagstat", "ds_1-HeLa-S3.bam.flagstat", "ds_2-HeLa-S3.bam.flagstat", "ds_GM12878-1-RNA.bam.flagstat", "dedup_GM12878-1-RNA-1.bam.flagstat", "ds_HeLs-S3.bam.flagstat" )
#rna_qc_all.sorted$flagstat.filename <- flagstat.filenames

#row.names(rna_qc_all_old_vs_new.sorted) <- rna_qc_all_old_vs_new.sorted$pretty.sample
row.names(rna_qc_all_2015_2014_2013.sorted) <-  rna_qc_all_2015_2014_2013.sorted$pretty.sample

old <- par()
par(mar=c(10,5,5,5))
#barplot (rna_qc_all_old_vs_new.sorted$rRNA.rate, col = c("lightgreen", "lightblue")[rna_qc_all_old_vs_new.sorted$center], names.arg = row.names(rna_qc_all_old_vs_new.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.004), ylab = "rRNA rate")
barplot (rna_qc_all_2015_2014_2013.sorted$rRNA.rate, col = c("lightblue", "blue", "darkblue")[rna_qc_all_2015_2014_2013.sorted$center], names.arg = row.names(rna_qc_all_2015_2014_2013.sorted), las = 2, cex.names  = 0.9, ylim = c(0,0.004), ylab = "rRNA rate")


# Strand specificity
rna_qc_all_old_vs_new.sorted$strand.spec <- rna_qc_all_old_vs_new.sorted$End.1.Sense / rna_qc_all_old_vs_new.sorted$End.2.Sense
rna_qc_all_old_vs_new.sorted$strand.spec.log10 <- log10(rna_qc_all_old_vs_new.sorted$strand.spec)
barplot(rna_qc_all_old_vs_new.sorted$strand.spec.log10, col = c("lightgreen", "lightblue")[rna_qc_all_old_vs_new.sorted$center], names.arg = row.names(rna_qc_all_old_vs_new.sorted), las = 2, cex.names  = 0.9, ylim = c(-2.5,0.1), ylab = "Strand specificity")

#PCR duplicate levels

## Correlations / scatterplots
old_new.cuff.GM12878 <- readCufflinks("../diff_out_GM12878_masked_separate_old_vs_new")
old_new.cuff.HeLa <- readCufflinks("../diff_out_HeLa_masked_separate_old_vs_new")

#GM12878
csScatterMatrix(genes(old_new.cuff.GM12878))
separate.fpkm.GM12878.old.new <- fpkmMatrix(genes(old_new.cuff.GM12878))
cor(separate.fpkm.GM12878.old.new, method = "spearman")

#HeLa-S3
csScatterMatrix(genes(old_new.cuff.HeLa))
separate.fpkm.HeLa.old.new <- fpkmMatrix(genes(old_new.cuff.HeLa))
cor(separate.fpkm.HeLa.old.new, method = "spearman")

