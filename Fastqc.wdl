## Task associated with fastqc for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task fastqc {

     String docker_image
     String sample_name
     File input_bam

     
     command <<<
      fastqc ~{input_bam}
     >>>
     runtime{
      docker: docker_image
     }
     output{
      File output_html = "${sample_name}_fastqc.html"
     }
}

task picard_metrics {
    
    String docker_image
    String sample_name
    File input_bam

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
     -I ~{input_bam} \
     -O ~{sample_name}-alignment.txt
     
     picard -Xmx{max_memory}G CollectInsertSizeMetrics \
     -I ~{input_bam} \
     -O ~{sample_name}-insert-size.txt \
     -H ~{sample_name}-hist.pdf 
     
    >>>
    runtime {
     docker: docker_image
    }
    output {
     File alignment   = "${sample_name}-alignment.txt"
     File insert_size = "${sample_name}-insert-size.txt"
     File hist        = "${sample_name}-hist.pdf"
    }
}


task qualimap {
    
    String docker_image
    File input_bam


    String memory
    String threads
    
    command {
      qualimap bamqc -bam ${input_bam} -nt ${threads} --collect-overlap-pairs --skip-duplicated --java-mem-size=${memory}G
    }
    runtime {
     docker: docker_image
    }
    output {
     File qualimap_report = "qualimapReport.html"
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

    
    command {
     MethylDackel mbias ${ref_fasta} ${bam_file}
    }
    runtime {
     docker: docker_image
    }
    output {
     File mbias_ot = "${sample_name}_OT.svg"
     File mbias_ob = "${sample_name}_OB.svg"
    }
    
}

task multiqc {

    String docker_image
    String sample_name
    File alnMetrics_input
    File insertMetrics_input
    File fastqc_input
    File qualimap_input

    
    command {
     echo ${alnMetrics_input} > my_file_list.txt
     echo ${insertMetrics_input} >> my_file_list.txt
     echo ${fastqc_input} >> my_file_list.txt
     echo ${qualimap_input} >> my_file_list.txt
     multiqc --file-list my_file_list.txt

    }
    runtime {
     docker: docker_image
    }
    output {
    File multiqc_report = "multiqc_report.html"
    }
    
}

task methylation_metrics {
    

    String docker_image
    File bed_graphs
    String wg_blimp_R_script_path
    
    command {
     R ${wg_blimp_R_script_path}/methylationMetrics.R
    }
    runtime {
     docker: docker_image
    }
    output {
     File methylation_metrics = "methylation_metrics.csv"
    }
    
}

