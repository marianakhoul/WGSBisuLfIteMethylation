## Task associated with alignment for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task bwameth_indexing {
	input{
		File ref_fasta
		String docker_image
		File bwameth_script
		File log
		String ref_fasta_name = basename(ref_fasta,".fa")

		String ref_amb = ref_fasta_name + ".amb"
	  	String ref_ann = ref_fasta_name + ".ann"
	  	String ref_bwt = ref_fasta_name + ".bwt"
	  	String ref_pac = ref_fasta_name + ".pac"
	  	String ref_sa = ref_fasta_name + ".sa"
	  	String ref_fasta_index = ref_fasta_name + ".fa.fai"
	}
	command {
		${bwameth_script} index ${ref_fasta} 2> ${log}
        samtools faidx ${ref_fasta} 2> ${log}
	}
	runtime {
		docker: docker_image
	}
	output{
		File ref_amb = ${ref_amb}
	  	File ref_ann = ${ref_ann}
	  	File ref_bwt = ${ref_bwt}
	  	File ref_pac = ${ref_pac}
	  	File ref_sa = ${ref_sa}
	  	File ref_fasta_index = ${ref_fasta_index}
	}
}
