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

task bwameth_align {

   File ref_amb
   File ref_ann
   File ref_bwt
   File ref_pac
   File ref_sa
   File ref_fasta_index
   File ref_fasta

   File fastq_file_1
   File fastq_file_2

   Int threads
   File bwameth_script
   String docker_image
   File log
   String alignment_dir
   String sample_name

   command {
    set -o pipefail
    set -e

    ${bwameth_script} -t ${threads} --reference ${ref_fasta} ${fastq_file_1} ${fastq_file_2} 2> ${log} \
    | \ 
    samtools view -b - > ${alignment_dir}${sample_name}.unsorted.bam 2>/dev/null
  }
  runtime {
    docker: docker_image
  }
  output{
    File output_unsorted_bam = "${alignment_dir}${sample_name}.unsorted.bam"
  }
}

task sort_bam {

   File input_bam
   Int threads
   String sample_name
   File log
   String alignment_dir
   String docker_image
  
   command <<<
    samtools sort -o ${alignment_dir}${sample_name}.sorted.bam -@ ${threads} ${input_bam} 2> {$log}
   >>>
   runtime {
    docker: docker_image
  }
   output{
    File output_sorted_bam = "${alignment_dir}${sample_name}.sorted.bam"
  }
}
