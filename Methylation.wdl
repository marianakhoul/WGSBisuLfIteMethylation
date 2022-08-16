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
import "https://raw.githubusercontent.com/marianakhoul/WGSBisuLfIteMethylation/main/Fastqc.wdl" as Fastqc
import "https://raw.githubusercontent.com/marianakhoul/WGSBisuLfIteMethylation/main/DMR_Calling.wdl" as DMR_Calling
import "https://raw.githubusercontent.com/marianakhoul/WGSBisuLfIteMethylation/main/DMR_Comparison.wdl" as DMR_Comparison
import "https://raw.githubusercontent.com/marianakhoul/WGSBisuLfIteMethylation/main/Segmentation.wdl" as Segmentation


workflow WGSBisuLfIteMethylation {
    
    # Parameters
    String wg_blimp_R_script_path
    String sample_name
    String case
    String control
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
    
    # Camel path
    String camel_modules_path

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
    call Alignment.bwameth_align {
        input:
            ref_fasta = ref_fasta,
            ref_amb = bwameth_indexing.ref_amb,
            ref_ann = bwameth_indexing.ref_ann,
            ref_bwt = bwameth_indexing.ref_bwt,
            ref_pac = bwameth_indexing.ref_pac,
            ref_sa = bwameth_indexing.ref_sa,
            ref_fasta_index = bwameth_indexing.ref_fasta_index,
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
    
    call Alignment.mark_duplicates {
        input:
            sample_name = sample_name,
            input_bam = sort_bam.output_sorted_bam,
            docker_image = picard_docker
    }
    call Alignment.index_bam {
        input:
            sample_name = sample_name,
            input_bam = mark_duplicates.output_bam,
            docker_image = gotc_docker
    }
    
    ## QC
    call Fastqc.fastqc {
        input:
            docker_image = fastqc_docker,
            sample_name = sample_name,
            input_bam = mark_duplicates.output_bam
    }
    
    call Fastqc.picard_metrics {
        input:
            ref_fasta = ref_fasta,
            ref_amb = bwameth_indexing.ref_amb,
            ref_ann = bwameth_indexing.ref_ann,
            ref_bwt = bwameth_indexing.ref_bwt,
            ref_pac = bwameth_indexing.ref_pac,
            ref_sa = bwameth_indexing.ref_sa,
            ref_fasta_index = bwameth_indexing.ref_fasta_index,
            input_bam = mark_duplicates.output_bam,
            sample_name = sample_name,
            docker_image = picard_docker
    }
    
    call Fastqc.qualimap {
        input:
            docker_image = qualimap_docker,
            sample_name = sample_name,
            input_bam = mark_duplicates.output_bam
    }
    
    call Fastqc.mbias {
        input:
            bam_index = index_bam.indexed_bam,
            bam_file = mark_duplicates.output_bam,
            ref_fasta = ref_fasta,
            ref_amb = bwameth_indexing.ref_amb,
            ref_ann = bwameth_indexing.ref_ann,
            ref_bwt = bwameth_indexing.ref_bwt,
            ref_pac = bwameth_indexing.ref_pac,
            ref_sa = bwameth_indexing.ref_sa,
            ref_fasta_index = bwameth_indexing.ref_fasta_index,
            sample_name = sample_name,
            docker_image = MethylDackel_docker
    }
    
    call Fastqc.multiqc {
        input:
            docker_image = multiqc_docker,
            alnMetrics_input = picard_metrics.alignment,
            insertMetrics_input = picard_metrics.insert_size,
            fastqc_input = fastqc.output_html,
            sample_name = sample_name,
            qualimap_input = qualimap.qualimap_report
    }
    
    call DMR_Calling.methyl_dackel {
        input:
            ref_fasta = ref_fasta,
            ref_amb = bwameth_indexing.ref_amb,
            ref_ann = bwameth_indexing.ref_ann,
            ref_bwt = bwameth_indexing.ref_bwt,
            ref_pac = bwameth_indexing.ref_pac,
            ref_sa = bwameth_indexing.ref_sa,
            ref_fasta_index = bwameth_indexing.ref_fasta_index,
            docker_image = MethylDackel_docker,
            bam_index = index_bam.indexed_bam,
            bam_file = mark_duplicates.output_bam,
            sample_name = sample_name
    }
    
    call Fastqc.methylation_metrics {
            input:
                    bed_graphs = methyl_dackel.methyl_dackel_output,
                    docker_image = R_docker
        }
    
    call DMR_Calling.bedgraph_to_methylation_ratio {
        input:
            methyl_dackel_output = methyl_dackel.methyl_dackel_output,
            sample_name = sample_name,
            docker_image = R_docker
    }
    
    call DMR_Calling.metilene_input {
        input:
            bedgraph_ratio = bedgraph_to_methylation_ratio.bedgraph_ratio,
            docker_image = bedtools_docker
    }
    
    call DMR_Calling.metilene {
        input:
            metilene_input_file = metilene_input.metilene_input_file,
            docker_image = metilene_docker
    }
    
    call DMR_Calling.bsseq {
        input:
            bedgraph_ratio = bedgraph_to_methylation_ratio.bedgraph_ratio,
            docker_image = R_docker,
            wg_blimp_R_script_path = wg_blimp_R_script_path
    }
    
    call DMR_Comparison.dmr_combination {
        input:
            wg_blimp_R_script_path = wg_blimp_R_script_path,
            docker_image = R_docker,
            ref_fasta_index = bwameth_indexing.ref_fasta_index
            
    }
    
    call DMR_Comparison.dmr_coverage {
        input:
            input_bam = mark_duplicates.output_bam,
            bam_bai = index_bam.indexed_bam,
            docker_image = mosdepth_docker,
            sample_name = sample_name,
            input_bed = dmr_combination.bed_output
    }
    
    call DMR_Comparison.dmr_annotation {
        input:
            biotypes = biotypes,
            wg_blimp_R_script_path = wg_blimp_R_script_path,
            docker_image = R_docker,
            coverages = dmr_coverage.regions_output,
            cgi_annotation_file = cgi_annotation_file,
            combined_dmrs = dmr_combination.csv_output,
            tss_distances = tss_distances,
            gene_annotation_file = gene_annotation_file,
            repeat_masker_annotation_file = repeat_masker_annotation_file
    }
    
    call Segmentation.methylseekr {
        input:
            wg_blimp_R_script_path = wg_blimp_R_script_path,
            docker_image = R_docker,
            cgi_annotation_file = cgi_annotation_file,
            gene_annotation_file = gene_annotation_file,
            repeat_masker_annotation_file = repeat_masker_annotation_file,
            sample_name = sample_name,
            tss_distances = tss_distances,
            methylation_table = methyl_dackel.methyl_dackel_output,
            ref_fasta = ref_fasta,
            fdr_cutoff = fdr_cutoff,
            methylation_cutoff = methylation_cutoff
            
    }
    
    
}
