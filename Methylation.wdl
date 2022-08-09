## This pipeline is the conversion of the Snakefile from https://github.com/MarWoes/wg-blimp.git
##
## This is pipeline #2 for the DNAnexus Pilot Study. Purpose of this script is for Whole Genome BisuLfIte sequencing ##Methylation analysis Pipeline.
##
##
## Workflow Description

import "Alignment.wdl" as Alignment

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

	File bwameth_script

	## ALIGNMENT
	call bwameth_indexing {
		input:
			docker_image = bwa_meth_docker,
			ref_fasta = ref_fasta,
			bwameth_script = bwameth_script
	}
}
