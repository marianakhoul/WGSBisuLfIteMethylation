version 1.1

import "./Alignment.wdl" as Alignment
import "./Fastqc.wdl" as Fastqc
import "./DMR_Calling.wdl" as DMR_Calling
import "./multiqc.wdl" as MultiQC
#import "./DMR_Comparison.wdl" as DMR_Comparison
#import "./Segmentation.wdl" as Segmentation

workflow WGSBisuLfIteMethylation {
    input {    
        # Parameters
        String sample_name
        Array[String] biotypes
        File cgi_annotation_file
        File repeat_masker_annotation_file
        File gene_annotation_file
        Array[Int] tss_distances
        Int fdr_cutoff
        Float methylation_cutoff

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
    }
    
    # Docker Images
    String MethylDackel_docker = "nfcore/methylseq"
    String python_docker = "python"
    String mosdepth_docker = "zlskidmore/mosdepth"
    String metilene_docker = "quay.io/biocontainers/metilene:0.2.8--h516909a_0"
    String R_docker = "dx://file-GG8Vvgj0J4yj9z6XFy19xPQv"
    String fastqc_docker = "pegi3s/fastqc"
    String qualimap_docker = "pegi3s/qualimap"  
    String bwa_meth_docker = "nfcore/methylseq"
    String gotc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
    String bedtools_docker = "biocontainers/bedtools"
    
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
            reference_fasta = reference_fasta,
            docker_image = bwa_meth_docker,
            bwameth_script = bwameth_script,
            sample_name = sample_name,
            fastq_file_1 = fastq_file_1,
            fastq_file_2 = fastq_file_2,
            threads = 8
    }
    
     call Alignment.sort_bam {
        input: 
            input_bam = bwameth_align.output_unsorted_bam,
            sample_name = sample_name,
            docker_image = gotc_docker,
            threads = 8
    }
    
    call Alignment.mark_duplicates {
        input:
            sample_name = sample_name,
            input_bam = sort_bam.output_sorted_bam,
            docker_image = gotc_docker
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
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            ref_fasta_index = ref_index,
            input_bam = mark_duplicates.output_bam,
            sample_name = sample_name,
            docker_image = gotc_docker
    }
    
    call Fastqc.qualimap {
        input:
            docker_image = qualimap_docker,
            input_bam = mark_duplicates.output_bam
    }
    
    call Fastqc.mbias {
        input:
            bam_index = index_bam.indexed_bam,
            bam_file = mark_duplicates.output_bam,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_index,
            ref_amb = ref_amb,
            ref_ann = ref_ann,
            ref_bwt = ref_bwt,
            ref_pac = ref_pac,
            ref_sa = ref_sa,
            sample_name = sample_name,
            docker_image = MethylDackel_docker
    }
    
    call MultiQC.multiqc as multiqc {
        input:
            input_logs = [
                picard_metrics.alignment,
                picard_metrics.insert_size,
                fastqc.output_html,
                qualimap.qualimap_report
                ]
    }
    
    call DMR_Calling.methyl_dackel {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_index,
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
    
    output {
        File output_unsorted_bam = bwameth_align.output_unsorted_bam
        File sorted_bam = sort_bam.output_sorted_bam
        File mark_duplicates_metrics = mark_duplicates.metrics
        File mark_duplicates_output_bam = mark_duplicates.output_bam
        File indexed_bam = index_bam.indexed_bam
        File picard_metrics_alignment   = picard_metrics.alignment
        File picard_metrics_insert_size = picard_metrics.insert_size
        File picard_metrics_hist        = picard_metrics.hist
        File fastqc_output_html = fastqc.output_html
        File qualimap_report = qualimap.qualimap_report
        File mbias_ot = mbias.mbias_ot
        File mbias_ob = mbias.mbias_ob
        File multiqc_report = multiqc.html_report
        File methyl_dackel_output = methyl_dackel.methyl_dackel_output
        File methylation_metrics_ouput = methylation_metrics.methylation_metrics_output
        File bedgraph_ratio = bedgraph_to_methylation_ratio.bedgraph_ratio
    } 
}
