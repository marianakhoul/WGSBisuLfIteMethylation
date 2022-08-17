## This pipeline is the conversion of the Snakefile from https://github.com/MarWoes/wg-blimp.git
##
## This is pipeline #2 for the DNAnexus Pilot Study. Purpose of this script is for Whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Workflow Description

#import "./Alignment.wdl" as Alignment
#import "./Fastqc.wdl" as Fastqc
#import "./DMR_Calling.wdl" as DMR_Calling
#import "./DMR_Comparison.wdl" as DMR_Comparison
#import "./Segmentation.wdl" as Segmentation

import "https://raw.githubusercontent.com/marianakhoul/WGSBisuLfIteMethylation/main/Alignment.wdl" as Alignment

workflow WGSBisuLfIteMethylation {
    
    # Parameters
    String wg_blimp_R_script_path
    String sample_name
    Array[String] biotypes
    File cgi_annotation_file
    File repeat_masker_annotation_file
    File gene_annotation_file
    Array[Int] tss_distances
    Int fdr_cutoff
    Float methylation_cutoff
    
    # Docker Images
    String MethylDackel_docker = "nfcore/methylseq"
    String python_docker = "python"
    String mosdepth_docker = "zlskidmore/mosdepth"
    String metilene_docker = "quay.io/biocontainers/metilene:0.2.8--h516909a_0"
    String R_docker = "r-base:latest"
    String multiqc_docker = "ewels/multiqc"
    String fastqc_docker = "pegi3s/fastqc"
    String picard_docker = "broadinstitute/picard"
    String qualimap_docker = "pegi3s/qualimap"  
    String bwa_meth_docker = "nfcore/methylseq"
    String gotc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    String bedtools_docker = "biocontainers/bedtools"
    
    # Reference Fasta
    File ref_fasta
    String reference_fasta

    # BWA Script
    File bwameth_script
    
    # Fastq files
    File fastq_file_1
    File fastq_file_2
    

    ## ALIGNMENT
    call Alignment.bwameth_indexing {
        input:
            docker_image = bwa_meth_docker,
            ref_fasta = ref_fasta,
            bwameth_script = bwameth_script
    }
}
