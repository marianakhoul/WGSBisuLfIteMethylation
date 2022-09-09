#!/usr/bin/env Rscript

library(optparse)
library(Rserve)
Rserve(args="--no-save")


option_list <- list(
                make_option(c("--bed_graphs"), type = "character", help = "bed_graphs file"),
                make_option(c("--output"), type = "character", help = "methylation_metrics output"),
                make_option(c("--meth_rate_on_chr"), type = "float",default="0.1", help = "methylation rate on chromosome"))

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

bed_graphs<-opt$bed_graphs
output<-opt$output
rate<-opt$meth_rate_on_chr


wgbs.annotateDMRs <- function (dmrFile, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, gzippedCoverageFiles, allowedBiotypes, promoterTSSDistances, annotatedDmrFile) {

  dmrs <- fread(dmrFile)

  dmrs <- annotation.annotateRegions(dmrs, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, allowedBiotypes, promoterTSSDistances)

  dmrLocationKeys <- paste0(dmrs$chr, ":", dmrs$start, "-", dmrs$end)

  coverages <- sapply(gzippedCoverageFiles, function (gzippedCoverageFile) {

    coverage = fread(cmd = paste("zcat", gzippedCoverageFile))
    coverageLocationKeys = paste0(coverage$V1, ":", coverage$V2, "-", coverage$V3)

    return(coverage$V4[match(dmrLocationKeys, coverageLocationKeys)])
  })

  dmrs$mean_cov <- apply(coverages, 1, mean)

  write.table(dmrs, file = annotatedDmrFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")

  return(dmrs)
}

wgbs.annotateDMRs(
     snakemake@input$combined_dmrs,
     snakemake@params$cgi_annotation_file,
     snakemake@params$gene_annotation_file,
     snakemake@params$repeat_masker_annotation_file,
     snakemake@input$coverages,
     snakemake@params$biotypes,
     snakemake@params$tss_distances,
     snakemake@output$annotated_dmrs
   )

