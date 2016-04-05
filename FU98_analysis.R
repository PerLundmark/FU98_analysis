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
rna_qc_all.sorted$raw.pcr.dup <- c(0.381183,0.394302,0.297067,0.281671,0.247557,0.217364,0.198524)
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

