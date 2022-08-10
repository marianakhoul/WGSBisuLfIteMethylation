## This pipeline is the conversion of the Snakefile from https://github.com/MarWoes/wg-blimp.git
##
## This is pipeline #2 for the DNAnexus Pilot Study. Purpose of this script is for Whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Workflow Description

import "/Users/marianakhoul/Desktop/DNAnexus_pipeline2/WGSBisuLfIteMethylation/Alignment.wdl" as Alignment
import "/Users/marianakhoul/Desktop/DNAnexus_pipeline2/WGSBisuLfIteMethylation/Fastqc.wdl" as Fastqc
import "/Users/marianakhoul/Desktop/DNAnexus_pipeline2/WGSBisuLfIteMethylation/DMR_Calling.wdl" as DMR_Calling
import "/Users/marianakhoul/Desktop/DNAnexus_pipeline2/WGSBisuLfIteMethylation/DMR_Comparison.wdl" as DMR_Comparison
import "/Users/marianakhoul/Desktop/DNAnexus_pipeline2/WGSBisuLfIteMethylation/Segmentation.wdl" as Segmentation


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
	
	# Docker Images
	String MethylDackel_docker
	String python_docker
	String R_docker = "r-base:latest"
	String multiqc_docker = "ewels/multiqc"
	String fastqc_docker = "pegi3s/fastqc"
	String picard_docker = "broadinstitute/picard"
	String qualimap_docker = "pegi3s/qualimap"
	String bwa_meth_docker = "pgcbioinfo/bwa-meth:latest"
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
	
	# Directories
	String alignment_dir
	String fastqc_dir
	String methylation_dir
	String mbias_dir
	String metilene_dir
	String camel_dir
	String bsseq_dir
	String dmr_dir

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
	  		ref_fasta = bwa_meth_docker,
	  		docker_image = python_docker,
	  		bwameth_script = bwameth_script,
	  		alignment_dir = alignment_dir,
	  		sample_name = sample_name,
	  		fastq_file_1 = fastq_file_1,
	  		fastq_file_2 = fastq_file_2
	}
	
	call Alignment.sort_bam {
		input: 
			input_bam = bwameth_align.output_unsorted_bam,
			sample_name = sample_name,
			alignment_dir = alignment_dir,
			docker_image = gotc_docker
	}
	
	call Alignment.mark_duplicates {
		input:
			sample_name = sample_name,
			input_bam = sort_bam.output_sorted_bam,
			alignment_dir = alignment_dir,
			docker_image = picard_docker
	}
	call Alignment.index_bam {
		input:
			sample_name = sample_name,
			input_bam = mark_duplicates.output_bam,
			alignment_dir = alignment_dir,
			docker_image = gotc_docker
	}
	
	## QC
	call Fastqc.fastqc {
		input:
			docker_image = fastqc_docker,
			sample_name = sample_name,
			alignment_dir = alignment_dir,
			fastqc_dir = fastqc_dir,
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
			alignment_dir = alignment_dir,
			fastqc_dir = fastqc_dir,
			sample_name = sample_name,
			docker_image = picard_docker
	}
	
	call Fastqc.qualimap {
		input:
			docker_image = qualimap_docker,
			sample_name = sample_name,
			fastqc_dir = fastqc_dir,
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
			mbias_dir = mbias_dir,
			sample_name = sample_name,
			docker_image = MethylDackel_docker
	}
	
	call Fastqc.multiqc {
		input:
			docker_image = multiqc_docker,
			alnMetrics_input = picard_metrics.alignment,
			insertMetrics_input = picard_metrics.insert_size,
			fastqc_input = fastqc.output_html,
			fastqc_dir = fastqc_dir,
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
			sample_name = sample_name,
			methylation_dir = methylation_dir
	}
	
	call Fastqc.methylation_metrics {
        	input:
            		bed_graphs = methyl_dackel.methyl_dackel_output,
			fastqc_dir = fastqc_dir,
			docker_image = R_docker
    	}
	
	call DMR_Calling.bedgraph_to_methylation_ratio {
		input:
			methylation_dir = methylation_dir,
			methyl_dackel_output = methyl_dackel.methyl_dackel_output,
			sample_name = sample_name,
			docker_image = R_docker
	}
	
	call DMR_Calling.metilene_input {
		input:
			bedgraph_ratio = bedgraph_to_methylation_ratio.bedgraph_ratio,
			docker_image = bedtools_docker,
			metilene_dir = metilene_dir,
			methylation_dir = methylation_dir
	}
	
	call DMR_Calling.metilene {
		input:
			metilene_dir = metilene_dir,
			metilene_input_file = metilene_input.metilene_input_file,
			docker_image = MethylDackel_docker
	}
	
	call DMR_Calling.camel_index {
		input:
			camel_modules_path = camel_modules_path,
			docker_image = python_docker,
			ref_fasta = ref_fasta,
			reference_fasta = reference_fasta
	}
	
	call DMR_Calling.camel_call {
		input:
			camel_modules_path = camel_modules_path,
			docker_image = python_docker,
			reference_output = camel_index.reference_output,
			input_bam = mark_duplicates.output_bam,
			bam_bai = index_bam.indexed_bam,
			sample_name = sample_name,
			camel_dir = camel_dir
	}
	
	call DMR_Calling.camel_dmr {
		input:
			camel_modules_path = camel_modules_path,
			docker_image = python_docker,
			reference_output = camel_index.reference_output,
			case = case,
			control = control,
			camel_dir = camel_dir
	}
	
	call DMR_Calling.bsseq {
		input:
			bedgraph_ratio = bedgraph_to_methylation_ratio.bedgraph_ratio,
			docker_image = R_docker,
			methylation_dir = methylation_dir,
			bsseq_dir = bsseq_dir,
			wg_blimp_R_script_path = wg_blimp_R_script_path
	}
	
	call DMR_Comparison.dmr_combination {
		input:
			wg_blimp_R_script_path = wg_blimp_R_script_path,
			docker_image = R_docker,
			dmr_dir = dmr_dir,
			ref_fasta_index = bwameth_indexing.ref_fasta_index
			
	}
	
	call DMR_Comparison.dmr_coverage {
		input:
			input_bam = mark_duplicates.output_bam,
			bam_bai = index_bam.indexed_bam,
			docker_image = MethylDackel_docker,
			dmr_dir = dmr_dir,
			sample_name = sample_name,
			input_bed = dmr_combination.bed_output
	}
	
	call DMR_Comparison.dmr_annotation {
		input:
			biotypes = biotypes,
			wg_blimp_R_script_path = wg_blimp_R_script_path,
			docker_image = R_docker,
			dmr_dir = dmr_dir,
			coverages = dmr_coverage.regions_output,
			cgi_annotation_file = cgi_annotation_file,
			combined_dmrs = dmr_combination.csv_output
	}
}
