## Task associated with alignment for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task bwameth_indexing {

    File ref_fasta
    String docker_image
    File bwameth_script
    File log
    String ref_fasta_name = basename(ref_fasta,".fa")

  
  command <<<
     ~{bwameth_script} index ~{ref_fasta} 2> ~{log}
    samtools faidx ~{ref_fasta} 2> ~{log}
  >>>
  runtime {
     docker: docker_image
  }
  output {
     File ref_amb = "${ref_fasta_name}.amb"
     File ref_ann = "${ref_fasta_name}.ann"
     File ref_bwt = "${ref_fasta_name}.bwt"
     File ref_pac = "${ref_fasta_name}.pac"
     File ref_sa = "${ref_fasta_name}.sa"
     File ref_fasta_index = "${ref_fasta_name}.fa.fai"
  }
}
