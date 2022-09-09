#!/usr/bin/env Rscript

library(optparse)
library(Rserve)
library(data.table)
Rserve(args="--no-save")


option_list <- list(
                make_option(c("--bed_graphs"), type = "character", help = "bed_graphs file"),
                make_option(c("--output"), type = "character", help = "methylation_metrics output"),
                make_option(c("--meth_rate_on_chr"), type = "integer",default="0.1", help = "methylation rate on chromosome"))

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

bed_graphs<-opt$bed_graphs
output<-opt$output
rate<-opt$meth_rate_on_chr


wgbs.methylation.computeMeanMethylation <- function(methylationTable) {

  methylatedValues <- methylationTable$V5
  allValues <- methylationTable$V5 + methylationTable$V6

  return(mean(methylatedValues / allValues))

}

wgbs.methylation.computeMethylationMetrics <- function(bedGraphFiles, methylationMetricsFile, conversionChromosomes) {

  methylationValues <- lapply(bedGraphFiles, fread, skip = 1)

  numTrailingCharacters <- nchar("_CpG.bedGraph")

  sampleNames <- basename(bedGraphFiles)
  sampleNames <- substr(sampleNames, 1, nchar(sampleNames) - numTrailingCharacters)

  names(methylationValues) <- sampleNames

  methylationMetrics <- data.table(
    sample = names(methylationValues),
    `methylation (total)` = sapply(methylationValues, wgbs.methylation.computeMeanMethylation)
  )

  if (length(conversionChromosomes) > 0) {
  
    conversionRates <- data.table(do.call(rbind, sapply(methylationValues, function (sampleValues) {
  
      chromosomeConversionRates <- sapply(conversionChromosomes, function(chr) {

        valuesOnChromosome <- sampleValues[V1 == chr]
  
        return(wgbs.methylation.computeMeanMethylation(valuesOnChromosome))
  
      })
  
      return(chromosomeConversionRates)
    }, simplify = FALSE)))

    colnames(conversionRates) <- paste("methylation (", conversionChromosomes, ")", sep = "")
    
    methylationMetrics <- cbind(methylationMetrics, conversionRates)
  }

  write.table(methylationMetrics, methylationMetricsFile, row.names = FALSE)

  return(methylationMetrics)
}


wgbs.methylation.computeMethylationMetrics(bed_graphs,output,rate)

