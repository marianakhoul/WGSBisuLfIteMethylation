#!/usr/bin/env Rscript

library(bsseq)
library(BiocParallel)
library(optparse)
library(Rserve)
Rserve(args="--no-save")


option_list <- list(
                make_option(c("--input"), type = "character", help = "bed_graphs file"),
                make_option(c("--rdata_file"), type = "character", help = "bed_graphs file"),
                make_option(c("--csv_file"), type = "character", help = "methylation_metrics output"),
                make_option(c("--pdf_file"), type = "float",default="0.1", help = "methylation rate on chromosome"))

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

rdata_file<-opt$rdata_file
csv_file<-opt$csv_file
pdf_file<-opt$pdf_file
input<-opt$input
threads<-10
min_cpg<-5
min_diff<-0.3

callDmrs <- function (methylDackelBedGraphFiles, sampleNames, group1Samples, group2Samples, threads, min_cpg, min_diff, localCorrect, rdatFile, csvFile, pdfFile) {

  biocParallelParam <- MulticoreParam(workers = threads)
  
  methylationData <- read.bismark(methylDackelBedGraphFiles,
                                  strandCollapse = FALSE,
                                  BPPARAM = biocParallelParam)
  
  sampleNames(methylationData) <- sampleNames

  # remove all NaN CpG loci: this should not happen in a WGBS experiment anyway, so the discarded
  # regions should be considered low-quality
  methylationData <- methylationData[!is.nan(rowSums(getMeth(methylationData, type = "raw")))]

  smoothedData <- BSmooth(methylationData, BPPARAM = biocParallelParam, verbose = TRUE)

  # add some color for later plotting
  pData <- pData(smoothedData)
  pData$col <- c(rep("red", length(group1Samples)), rep("blue", length(group2Samples)))
  pData(smoothedData) <- pData

  # filter out rows containing NA's
  invalidRows <- is.na(rowSums(getMeth(smoothedData)))

  print(paste("Filtering out", sum(invalidRows), "rows containing NA"))

  tstats <- BSmooth.tstat(smoothedData[!invalidRows],
                          group1 = group1Samples,
                          group2 = group2Samples,
                          estimate.var = "same",
                          local.correct = localCorrect,
                          mc.cores = threads)

  statType <- ifelse(localCorrect, "tstat.corrected", "tstat")

  dmrs0 <- dmrFinder(tstats, stat = statType)

  write.table(dmrs0, file = csvFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

  pdf(file = pdfFile, width = 10, height = 5)

  dmrs <- subset(dmrs0, n >= min_cpg & abs(meanDiff) >= min_diff)

  maxNumOfRegions <- ifelse(nrow(dmrs) > 100, 100, nrow(dmrs))
  plotManyRegions(smoothedData[!invalidRows], regions = dmrs[seq_len(maxNumOfRegions),], addRegions = dmrs, extend = 5000)

  dev.off()
  save.image(file = rdatFile)
}

callDmrs(input,
           snakemake@config$samples,
           snakemake@config$group1,
           snakemake@config$group2,
           threads,
           min_cpg,
           min_diff,
           snakemake@params$local_correct,
           rdata_file,
           csv_file,
           pdf_file)



