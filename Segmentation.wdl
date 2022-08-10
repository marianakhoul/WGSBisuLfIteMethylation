## Task associated with Segmentation for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task methylseekr {

    String wg_blimp_R_script_path
    string docker_image
    String sample_name
    File cgi_annotation_file
    File gene_annotation_file
    File repeat_masker_annotation_file
    String segmentation_dir
    Array[String] biotypes
    Array[Int] tss_distances
    Int min_cov
    File methylation_table
    File ref_fasta
    
    command {
      R ${wg_blimp_R_script_path}/methylseekRSegmentation.R
    }
    runtime {
     docker: docker_image
    }
    output {
     File pmd_all = "${segmentation_dir}/pmd-all.csv"
     File umr_lmr_all = "${segmentation_dir}/umr-lmr-all.csv"
     File pmd_segments = "${segmentation_dir}/${sample_name}/pmd-segments.csv"
     File umr_lmr_with_pmd = "${segmentation_dir}/${sample_name}/LMRUMRwithPMD/umr-lmr.csv"
     File umr_lmr_without_pmd = "${segmentation_dir}/${sample_name}/LMRUMRwithoutPMD/umr-lmr.csv"
    }
}
