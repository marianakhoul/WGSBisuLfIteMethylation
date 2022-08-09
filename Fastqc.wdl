## Task associated with fastqc for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task fastqc {

     String docker_image
     String sample_name
     File input_bam
     String alignment_dir
     String fastqc_dir
     File log
     
     command <<<
      fastqc -o ~{fastqc_dir} ~{alignment_dir}~{input_bam} &> ~{log}
     >>>
     runtime{
      docker: docker_image
     }
     output{
      File output_html = "${fastqc_dir}/fastqc/${sample_name}_fastqc.html"
     }
}
