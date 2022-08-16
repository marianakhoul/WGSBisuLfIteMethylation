## Task associated with DMR comparison for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task dmr_combination {
    

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
     File bed_output = "combined-dmrs.csv"
     File csv_output = "combined-dmrs.bed"
    }
}


task dmr_coverage {
    

    File input_bam
    File bam_bai
    File input_bed
    String sample_name
    Int annotation_min_mapq
    Int threads
    String docker_image
    
    command {
     mosdepth --threads ${threads} --no-per-base --mapq ${annotation_min_mapq} --by ${input_bed} ${input_bam}

    }
    runtime {
     docker: docker_image
    }
    output {
     File regions_output = "${sample_name}.regions.bed.gz"
    }
}


task dmr_annotation {
    
    String wg_blimp_R_script_path
    String docker_image

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
     File annotated_dmrs = "annotated-dmrs.csv"
    }
}
