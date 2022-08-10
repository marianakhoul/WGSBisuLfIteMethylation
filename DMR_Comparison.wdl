## Task associated with DMR comparison for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task dmr_combination {
    
    String dmr_dir
    String docker_image
    String wg_blimp_R_script_path
    File ref_fasta_index
    
    command {
     R ${wg_blimp_R_script_path}/dmrCombination.R
    }
    runtime {
     docker: docker_image
    }
    output {
     File bed_output = "${dmr_dir}/combined-dmrs.csv"
     File csv_output = "${dmr_dir}/dmr-coverage/combined-dmrs.bed"
    }
}


task dmr_coverage {
    
    String dmr_dir
    File input_bam
    File bam_bai
    File input_bed
    String sample_name
    Int annotation_min_mapq
    Int threads
    String docker_image
    
    command {
     mosdepth --threads ${threads} --no-per-base --mapq ${annotation_min_mapq} --by ${input_bed} ${dmr_dir}/${sample_name} ${input_bam}

    }
    runtime {
     docker: docker_image
    }
    output {
     File regions_output = "${dmr_dir}/dmr-coverage/${sample_name}.regions.bed.gz"
    }
}


task dmr_annotation {
    
    String wg_blimp_R_script_path
    String docker_image
    String dmr_dir
    File coverages
    File combined_dmrs
    Array[String] biotypes
    Array[Int] tss_distances
    File cgi_annotation_file
    File gene_annotation_file
    File repeat_masker_annotation_file
    
    command {
     R ${wg_blimp_R_script_path}/dmrAnnotation.R
    }
    runtime {
     docker: docker_image
    }
    output {
     File annotated_dmrs = "${dmr_dir}/annotated-dmrs.csv"
    }
}
