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
    File ref_amb
 	File ref_ann
 	File ref_bwt
 	File ref_pac
 	File ref_sa
 	File ref_index
    File reference_fasta

    # BWA Script
    File bwameth_script
    
    # Fastq files
    File fastq_file_1
    File fastq_file_2
    

    ## ALIGNMENT
    call Alignment.bwameth_align {
        input:
            ref_fasta = ref_fasta,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            ref_fasta_index = ref_index,
            reference_fasta=reference_fasta
            docker_image = bwa_meth_docker,
            bwameth_script = bwameth_script,
            sample_name = sample_name,
            fastq_file_1 = fastq_file_1,
            fastq_file_2 = fastq_file_2
    }
    
     call Alignment.sort_bam {
        input: 
            input_bam = bwameth_align.output_unsorted_bam,
            sample_name = sample_name,
            docker_image = gotc_docker
    }
    
    output {
        File output_unsorted_bam = bwameth_align.output_unsorted_bam
        File sorted_bam = sort_bam.output_sorted_bam
  }
}
