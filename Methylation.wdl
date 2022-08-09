version 1.0

## This pipeline is the conversion of the Snakefile from https://github.com/MarWoes/wg-blimp.git
##
## This is pipeline #2 for the DNAnexus Pilot Study. Purpose of this script is for Whole Genome BisuLfIte sequencing ##Methylation analysis Pipeline.
##
##
## Workflow Description

import ./Alignment.wdl as Alignment

workflow WGSBisuLfIteMethylation {
	
	String wg_blimp_R_script_path = "./scripts"
	String sample_name

	String python_docker
	
	# Reference Fasta
	File ref_fasta

	File bwameth_script

	## ALIGNMENT
	call bwameth_indexing {
		input:
			docker_image = python_docker,
			ref_fasta = ref_fasta,
			bwameth_script = bwameth_script
	}
}
