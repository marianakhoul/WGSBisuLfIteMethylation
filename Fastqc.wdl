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
      fastqc -o ~{fastqc_dir}/fastqc ~{alignment_dir}~{input_bam} &> ~{log}
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


task qualimap {
    
    String docker_image
    String sample_name
    File input_bam
    String fastqc_dir
    File log
    String memory
    String threads
    
    command {
      qualimap bamqc -bam ${input_bam} -outdir ${fastqc_dir}/qualimap -nt ${threads} --collect-overlap-pairs --skip-duplicated --java-mem-size=${memory}G &> ${log}
    }
    runtime {
     docker: docker_image
    }
    output {
     File qualimap_report = "${fastqc_dir}/qualimap/${sample_name}/qualimapReport.html"
    }
}
task mbias {

    String docker_image
    String sample_name
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fasta_index
    File ref_fasta
    File bam_index
    File bam_file
    String mbias_dir
    File log
    
    command {
     MethylDackel mbias ${ref_fasta} ${bam_file} ${mbias_dir}/${sample_name} &> ${log}
    }
    runtime {
     docker: docker_image
    }
    output {
     File mbias_ot = "${mbias_dir}/${sample_name}_OT.svg"
     File mbias_ob = "${mbias_dir}/${sample_name}_OB.svg"
    }
    
}

task multiqc {

    String docker_image
    String sample_name
    File alnMetrics_input
    File insertMetrics_input
    File fastqc_input
    File qualimap_input
    String fastqc_dir
    File log
    
    command {
     multiqc -f -o ${fastqc_dir} ${fastqc_dir}/fastqc ${fastqc_dir}/picard-metrics ${fastqc_dir}/qualimap &> ${log} 

    }
    runtime {
     docker: docker_image
    }
    output {
    File multiqc_report = "${fastqc_dir}/multiqc_report.html"
    }
    
}

task methylation_metrics {
    
    String fastqc_dir
    File log
    File bed_graphs
    
    command {
     scripts/methylationMetrics.R
    }
    runtime {
    }
    output {
     File methylation_metrics = "${fastqc_dir}/methylation_metrics.csv"
    }
    
}

