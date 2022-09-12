#!/usr/bin/env Rscript

library(optparse)
library(Rserve)
library(data.table)
Rserve(args="--no-save")


option_list <- list(
                make_option(c("--input"), type = "character", help = "MethylDackel bed graph output file as input"),
                make_option(c("--output"), type = "character", help = "bedgraph_to_methylation_ratio {sample_name}_CpG_ratio.bedGraph output file"))

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

input<-opt$input
output<-opt$output


wgbs.metilene.transformMethylationValues <- function (inputFile, outputFile) {

  fileContent <- fread(inputFile, data.table = FALSE, stringsAsFactors = FALSE)

  # see https://github.com/dpryan79/MethylDackel#single-cytosine-methylation-metrics-extraction
  methylatedReadsColumnIndex   <- 5
  unmethylatedReadsColumnIndex <- 6

  methylationRatios <- fileContent[,methylatedReadsColumnIndex] / (fileContent[,methylatedReadsColumnIndex] + fileContent[,unmethylatedReadsColumnIndex])

  transformedValues <- data.frame(
    chr = fileContent[,1],
    start = fileContent[,2],
    end = fileContent[,3],
    value = methylationRatios,
    stringsAsFactors = FALSE
  )

  write.table(transformedValues, file = outputFile, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
}


wgbs.metilene.transformMethylationValues(input,output)

