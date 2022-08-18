## Task associated with alignment for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

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

   String sample_name
   String docker_image

   command {
    set -o pipefail
    set -e

    ${bwameth_script} -t ${threads} --reference ${ref_fasta} ${fastq_file_1} ${fastq_file_2} \
    | \ 
    samtools view -b - > ${sample_name}.unsorted.bam 2>/dev/null
  }
  runtime {
    docker: docker_image
  }
  output{
    File output_unsorted_bam = "${sample_name}.unsorted.bam"
  }
}

task sort_bam {

   File input_bam
   Int threads
   String sample_name

   String docker_image
  
   command <<<
    samtools sort -o ${sample_name}.sorted.bam -@ ${threads} ${input_bam}
   >>>
   runtime {
    docker: docker_image
  }
   output{
    File output_sorted_bam = "${sample_name}.sorted.bam"
  }
}

task mark_duplicates {

      String sample_name
      File input_bam
      Float max_memory
      String docker_image

     command {
      picard -Xmx${max_memory}G \
      MarkDuplicates \
      -I ${input_bam} \
      -O ${sample_name}.bam \
      -M ${sample_name}-dup-metrics.txt 
     }
     runtime {
      docker: docker_image
     }
     output{
      File output_bam = "${sample_name}.bam"
      File metrics = "${sample_name}-dup-metrics.txt" 
     }
}

task index_bam {

     String sample_name
     File input_bam
  
     String docker_image

     command {
		samtools index ${input_bam} ${sample_name}.bai
     }
     runtime {
		docker: docker_image
     }
     output{
		File indexed_bam = "${sample_name}.bai"
     }
}
