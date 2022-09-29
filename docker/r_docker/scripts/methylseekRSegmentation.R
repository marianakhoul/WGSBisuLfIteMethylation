#!/usr/bin/env Rscript

library(optparse)
library(MethylSeekR)
library(tidyverse)
library(data.table)
library(rtracklayer)
library(parallel)
library(Biostrings)


option_list <- list(
                make_option(c("--input_ref"), type = "character", help = "input reference fasta"),
                make_option(c("--cgiAnnotation"), type = "character", help = "cgi annotation file path"),
                make_option(c("--geneAnnotation"), type = "character", help = "gene annotation file path"),
                make_option(c("--repeatMaskerAnnotation"), type = "character", help = "repeat masker annotation file path"),
                make_option(c("--threads"), type = "int", default = "10",help = "number of threads to use"),
                make_option(c("--tss_distances"), type = "character", help = "tss distance array"),
                make_option(c("--methylation_cutoff"), type = "float",default = "0.5", help = "methylation cutoff number"),
                make_option(c("--fdr_cutoff"), type = "int",default = "5", help = "FDR cutoff number"),
                make_option(c("--biotypes"), type = "character", help = "allowed biotypes array"),
                make_option(c("--pmd_all"), type = "character", help = "pmd_all output file name"),
                make_option(c("--umr_lmr_all"), type = "character", help = "umr_lmr_all output file name"),
                make_option(c("--min_coverage"), type = "int",default = "5", help = "Minimum coverage number"),
                


)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)


fastaRefFile <- opt$input_ref
cgiAnnotationFile <- opt$cgiAnnotation
geneAnnotationFile <- opt$geneAnnotation
repeatMaskerAnnotationFile <- opt$repeatMaskerAnnotation
samples <- snakemake@params$samples
targetDir <- "./"
methylationDir <- "./"
calibrationChr <- snakemake@params$calibration_chr
allowedBiotypes <- opt$biotypes
tssDistances <- opt$tss_distances
numThreads <- opt$threads
targetPmdFile <- opt$pmd_all
targetUmrLmrFile <- opt$umr_lmr_all
minimumCoverage <- opt$min_coverage
fdrCutoff <- opt$fdr_cutoff
methylationCutoff <- opt$methylation_cutoff

### GLOBALS

SEGMENTATION_WIDTH_PX <- 1366
SEGMENTATION_HEIGHT_PX <- 768
SEGMENTATION_FONT_SIZE <- 14

### FUNCTIONS
segmentIntoUMRsAndLMRs <- function(sample, methylationRanges, cgiRanges, numThreads, targetDir, genomeSeq, seqLengths, pmdSegments, sequencesStartingWithChr, minCoverage, fdrCutoff, methylationCutoff) {

  png(paste0(targetDir, "fdr.stats.png"))
  fdrStats <- calculateFDRs(
    m = methylationRanges,
    CGIs = cgiRanges,
    PMDs = pmdSegments,
    num.cores = numThreads,
    minCover = minCoverage
  )
  dev.off()

  selectedN <- as.integer(names(fdrStats$FDRs[as.character(methylationCutoff), ][fdrStats$FDRs[as.character(methylationCutoff), ] < fdrCutoff])[1])

  png(paste0(targetDir, "umr-lmr-heatmap.png"))
  segments <- segmentUMRsLMRs(
    m = methylationRanges,
    meth.cutoff = methylationCutoff,
    nCpG.cutoff = selectedN,
    PMDs = pmdSegments,
    num.cores = numThreads,
    myGenomeSeq = genomeSeq,
    seqLengths = seqLengths,
    minCover = minCoverage
  )
  dev.off()

  png(paste0(targetDir, "umr-lmr-segmentation.png"), width = SEGMENTATION_WIDTH_PX, height = SEGMENTATION_HEIGHT_PX, pointsize = SEGMENTATION_FONT_SIZE)
  plotFinalSegmentation(
    m = methylationRanges,
    segs = segments,
    PMDs = pmdSegments,
    meth.cutoff = methylationCutoff,
    minCover = minCoverage
  )
  dev.off()

  umrLmrTable <- data.table(
    sample = sample,
    chr = as.character(seqnames(segments)),
    start = start(segments),
    end = end(segments),
    type = values(segments)$type,
    num_cpgs_filtered = values(segments)$nCG.segmentation,
    num_cpgs_ref = values(segments)$nCG,
    mean_methylation = values(segments)$pmeth,
    median_methylation = values(segments)$median.meth,
    pmds_masked = !(length(pmdSegments) == 1 && is.na(pmdSegments))
  )

  fwrite(umrLmrTable, file = paste0(targetDir, "umr-lmr.csv"))
}

# HACK: create own class for DNAStringSets with MethylSeekR. Otherwise one would always have to
# install own packages. Also, general FASTA files could not be used otherwise
setClass("BlimpDnaStringSet", contains = "DNAStringSet")
setMethod("getSeq", "BlimpDnaStringSet", function (x, names, ...) {
  callNextMethod(x, names)
})



snakemake@source("regionAnnotation.R")

### SCRIPT

fastaRef <- readDNAStringSet(fastaRefFile)

# Remove description of sequence names according to FASTA format specification.
# Achieved by cutting everything off that comes after the first whitespace
names(fastaRef) <- sapply(str_split(names(fastaRef), "\\s"), `[`, 1)

class(fastaRef) <- "BlimpDnaStringSet"

if (is.null(cgiAnnotationFile) || !file.exists(cgiAnnotationFile)) {

  stop("[ERROR] No CGIs set for segmentation. CGIs must be set to perform segmentation with MethylSeekR.")

}

cgiAnnotation <- fread(cmd = paste("zcat", cgiAnnotationFile))

# detect chromosome naming ('chr1' vs '1') and
# adjust cgi chromosome names accordingly
if(sum(str_count(names(fastaRef), "^chr")) == 0) {
  cgiAnnotation$chrom <- str_remove(cgiAnnotation$chrom, "^chr")
  calibrationChr <- str_remove(calibrationChr, "^chr")
}

cgiRanges <- GRanges(cgiAnnotation$chrom, IRanges(cgiAnnotation$chromStart, cgiAnnotation$chromEnd))
cgiRanges <- suppressWarnings(resize(cgiRanges, 5000, fix = "center"))



for (sample in samples) {

  methylationValues <- fread(
    paste0(methylationDir, sample, "_CpG.bedGraph"),
    skip = 1L,
    header = FALSE,
    col.names = c("chr", "start", "end", "methPercent", "methylated", "unmethylated")
  )

  sampleRanges <- GRanges(methylationValues$chr, ranges = methylationValues$start)
  values(sampleRanges)$T <- methylationValues$methylated + methylationValues$unmethylated
  values(sampleRanges)$M <- methylationValues$methylated

  png(paste0(targetDir, sample, "/alphaCalibration.png"))
  pmdSegments <- segmentPMDs(
    m = sampleRanges,
    chr.sel = calibrationChr,
    num.cores = numThreads,
    seqLengths = seqlengths(fastaRef)
  )
  dev.off()

  png(paste0(targetDir, sample, "/alphaPMDRemoved.png"))
  plotAlphaDistributionOneChr(
    m = subsetByOverlaps(sampleRanges, pmdSegments[values(pmdSegments)$type == "notPMD"]),
    chr.sel = calibrationChr,
    num.cores = numThreads
  )
  dev.off()

  png(paste0(targetDir, sample, "/pmd.png"), width = SEGMENTATION_WIDTH_PX, height = SEGMENTATION_HEIGHT_PX, pointsize = SEGMENTATION_FONT_SIZE)
  plotPMDSegmentation(
    m = sampleRanges,
    segs = pmdSegments,
    minCover = minimumCoverage
  )
  dev.off()

  pmdSegmentTable <- data.table(
    sample = sample,
    chr = as.character(seqnames(pmdSegments)),
    start = start(pmdSegments),
    end = end(pmdSegments),
    type = values(pmdSegments)$type,
    num_cpgs = values(pmdSegments)$nCG
  )
  pmdSegmentTable <- pmdSegmentTable[type == "PMD"]

  fwrite(pmdSegmentTable, file = paste0(targetDir, sample, "/pmd-segments.csv"))

  segmentIntoUMRsAndLMRs(
    sample = sample,
    methylationRanges = sampleRanges,
    cgiRanges = cgiRanges,
    numThreads = numThreads,
    targetDir = paste0(targetDir, sample, "/LMRUMRwithPMD/"),
    genomeSeq = fastaRef,
    seqLengths = seqlengths(fastaRef),
    pmdSegments = pmdSegments,
    minCoverage = minimumCoverage,
    fdrCutoff = fdrCutoff,
    methylationCutoff = methylationCutoff
  )

  segmentIntoUMRsAndLMRs(
    sample = sample,
    methylationRanges = sampleRanges,
    cgiRanges = cgiRanges,
    numThreads = numThreads,
    targetDir = paste0(targetDir, sample, "/LMRUMRwithoutPMD/"),
    genomeSeq = fastaRef,
    seqLengths = seqlengths(fastaRef),
    pmdSegments = NA,
    minCoverage = minimumCoverage,
    fdrCutoff = fdrCutoff,
    methylationCutoff = methylationCutoff
  )

}

pmdFiles <- list.files(targetDir, pattern = "pmd-segments.csv", recursive = TRUE, full.names = TRUE)
umrLmrFiles <- list.files(targetDir, pattern = "umr-lmr.csv", recursive = TRUE, full.names = TRUE)

pmdContents <- lapply(pmdFiles, fread)
umrLmrContents <- lapply(umrLmrFiles, fread)

combinedPmdTable <- rbindlist(pmdContents)
combinedUmrLmrTable <- rbindlist(umrLmrContents)

annotatedPmdTable <- annotation.annotateRegions(
  combinedPmdTable,
  cgiAnnotationFile,
  geneAnnotationFile,
  repeatMaskerAnnotationFile,
  allowedBiotypes,
  tssDistances
)

annotatedUmrLmrTable <- annotation.annotateRegions(
  combinedUmrLmrTable,
  cgiAnnotationFile,
  geneAnnotationFile,
  repeatMaskerAnnotationFile,
  allowedBiotypes,
  tssDistances
)

fwrite(annotatedPmdTable, targetPmdFile)
fwrite(annotatedUmrLmrTable, targetUmrLmrFile)
