## Task associated with DMR comparison for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task dmr_combination {
    
    String docker_image
    File ref_index
    File metilene_input_file
    File bsseq_input

    command {
     Rscript /usr/local/bin/dmrCombination.R --bsseq_input ${bsseq_input} --metilene_input ${metilene_input_file} --fasta_index ${ref_index} --bed_output "combined-dmrs.csv" --csv_output 
"combined-dmrs.bed"
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
    
    String docker_image
    File coverages
    File combined_dmrs
    Array[String] biotypes
    Array[Int] tss_distances
    File cgi_annotation_file
    File gene_annotation_file
    File repeat_masker_annotation_file
    
    command {
     Rscript /usr/local/bin/dmrAnnotation.R
    }
    runtime {
     docker: docker_image
    }
    output {
     File annotated_dmrs = "annotated-dmrs.csv"
    }
}
