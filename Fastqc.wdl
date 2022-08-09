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

task picard_metrics {
    
    String docker_image
    String sample_name
    File input_bam
    String alignment_dir
    String fastqc_dir
    File log
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fasta_index
    File ref_fasta
    Int max_memory
    
    command <<<
    
     picard -Xmx{max_memory}G CollectAlignmentSummaryMetrics \
     -R ~{ref_fasta} \
     -I ~{alignment_dir}~{input_bam} \
     -O ~{fastqc_dir}/picard-metrics/~{sample_name}-alignment.txt &> ~{log}
     
     picard -Xmx{max_memory}G CollectInsertSizeMetrics \
     -I ~{alignment_dir}~{input_bam} \
     -O ~{fastqc_dir}/picard-metrics/~{sample_name}-insert-size.txt \
     -H ~{fastqc_dir}/picard-metrics/~{sample_name}-hist.pdf  &> ~{log}
     
    >>>
    runtime {
     docker: docker_image
    }
    output {
     File alignment   = "${fastqc_dir}/picard-metrics/${sample_name}-alignment.txt"
     File insert_size = "${fastqc_dir}/picard-metrics/${sample_name}-insert-size.txt"
     File hist        = "${fastqc_dir}/picard-metrics/${sample_name}-hist.pdf"
    }
}
