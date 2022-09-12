#!/usr/bin/env Rscript

library(GenomicRanges)
library(data.table)
library(tidyverse)
library(optparse)
library(Rserve)
Rserve(args="--no-save")


option_list <- list(
                make_option(c("--bsseq_input"), type = "character", help = "bsseq output as input file"),
                make_option(c("--bed_output"), type = "character", help = "bed file output"),
                make_option(c("--csv_output"), type = "character", help = "csv output file"),
                make_option(c("--csv_output"), type = "character", help = "csv output file"),)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

bsseq_input<-opt$bsseq_input
bed_output<-opt$bed_output
csv_output<-opt$csv_output
min_cpg<-5
min_diff<-0.3


wgbs.rangesFromData <- function (chr, start, end, num_cpgs, diff, tool, qValues) {


  toolRanges <- GRanges(chr, ranges = IRanges(start, end))
  toolRanges$num_cpgs <- num_cpgs
  toolRanges$tool <- tool
  toolRanges$diff <- diff
  toolRanges$q <- qValues

  return(toolRanges)

}

wgbs.loadBsseqData <- function (bsseqFile, minCpG, minDiff) {

  if (is.null(bsseqFile)) return(GRanges())

  rawData <- fread(bsseqFile, stringsAsFactors = FALSE, data.table = FALSE)

  validRows <- abs(rawData$meanDiff) >= minDiff & rawData$n >= minCpG

  validData <- rawData[validRows,]

  return(wgbs.rangesFromData(
    validData$chr,
    validData$start,
    validData$end,
    validData$n,
    validData$meanDiff,
    "bsseq",
    NA)
  )

}

wgbs.loadCamelData <- function (camelFile) {

  if (is.null(camelFile)) return(GRanges())

  rawData <- fread(camelFile, stringsAsFactors = FALSE, data.table = FALSE)

  colnames(rawData) <- c("chr", "start", "stop", "n", "core_diff")

  return(wgbs.rangesFromData(
    rawData$chr,
    rawData$start,
    rawData$stop,
    rawData$n,
    rawData$core_diff,
    "camel",
    NA)
  )
}

wgbs.loadMetileneData <- function (metileneFile) {

  if (is.null(metileneFile)) return(GRanges())

  rawData <- fread(metileneFile, stringsAsFactors = FALSE, data.table = FALSE)

  colnames(rawData) <- c("chr", "start", "stop", "q-value", "mean_diff", "num cpgs", "p (MWU)", "p (2D KS)", "g1 mean", "g2 mean")

  return(wgbs.rangesFromData(
    rawData$chr,
    rawData$start,
    rawData$stop,
    rawData$`num cpgs`,
    rawData$mean_diff,
    "metilene",
    rawData$`q-value`)
  )

}

wgbs.combineDmrs <- function (bsseqFile, camelFile, metileneFile, fastaIndex, csvOutput, bedOutput, minCpG, minDiff) {

  bsseqRanges    <- wgbs.loadBsseqData(bsseqFile, minCpG, minDiff)
  camelRanges    <- wgbs.loadCamelData(camelFile)
  metileneRanges <- wgbs.loadMetileneData(metileneFile)

  combinedRanges <- c(bsseqRanges, camelRanges, metileneRanges)

  seqs <- as.character(seqnames(combinedRanges))

  combinedDMRs <- data.frame(
    chr = seqs,
    start = start(combinedRanges),
    end = end(combinedRanges),
    num_cpgs = combinedRanges$num_cpgs,
    diff = combinedRanges$diff,
    qValue = combinedRanges$q,
    tool = combinedRanges$tool,
    stringsAsFactors = FALSE
  )

  fastaIndex <- fread(fastaIndex)

  # define own sort order
  seqFactors <- factor(seqs, levels = fastaIndex$V1)

  sortedDMRs <- combinedDMRs[order(seqFactors, combinedDMRs$start, combinedDMRs$end),]

  write.table(sortedDMRs, file = csvOutput, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")
  write.table(sortedDMRs[,c("chr", "start", "end")], file = bedOutput, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

  return(sortedDMRs)
}

 wgbs.combineDmrs(
    bsseq_input,
    #snakemake@input$camel,
    snakemake@input$metilene,
    snakemake@input$fasta_index,
    csv_output,
    bed_output,
    min_cpg,
    min_diff
  )

