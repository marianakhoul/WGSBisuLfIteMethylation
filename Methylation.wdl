## This pipeline is the conversion of the Snakefile from https://github.com/MarWoes/wg-blimp.git
##
## This is pipeline #2 for the DNAnexus Pilot Study. Purpose of this script is for Whole Genome BisuLfIte sequencing ##Methylation analysis Pipeline.
##
##
## Workflow Description

import "/Users/marianakhoul/Desktop/DNAnexus_pipeline2/WGSBisuLfIteMethylation/Alignment.wdl" as Alignment


workflow WGSBisuLfIteMethylation {
	
	String wg_blimp_R_script_path = "./scripts"
	String sample_name

	String multiqc_docker = "ewels/multiqc"
	String fastqc_docker = "pegi3s/fastqc"
	String picard_docker = "broadinstitute/picard"
	String qualimap_docker = "pegi3s/qualimap"
	String bwa_meth_docker = "pgcbioinfo/bwa-meth:latest"
	
	# Reference Fasta
	File ref_fasta

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
			ref_amb = Alignment.bwameth_indexing.ref_amb,
			ref_ann = Alignment.bwameth_indexing.ref_ann,
	  		ref_bwt = Alignment.bwameth_indexing.ref_bwt,
	  		ref_pac = Alignment.bwameth_indexing.ref_pac,
	  		ref_sa = Alignment.bwameth_indexing.ref_sa,
	  		ref_fasta_index = Alignment.bwameth_indexing.ref_fasta_index,
	  		ref_fasta = bwa_meth_docker,
	  		docker_image = python_docker,
	  		bwameth_script = bwameth_script,
	  		alignment_dir = alignment_dir,
	  		sample_name = sample_name,
	  		fastq_file_1 = fastq_file_1,
	  		fastq_file_2 = fastq_file_2
	}
}
