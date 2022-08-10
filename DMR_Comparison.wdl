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
