## Task associated with Segmentation for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task methylseekr {

    String docker_image
    String sample_name
    File cgi_annotation_file
    File gene_annotation_file
    File repeat_masker_annotation_file
    Array[String] biotypes
    Array[Int] tss_distances
    Int min_cov
    File methylation_table
    File ref_fasta
    Int fdr_cutoff
    Float methylation_cutoff
    
    command {
      Rscript /usr/local/bin/methylseekRSegmentation.R --cgiAnnotation cgi_annotation_file --geneAnnotation gene_annotation_file --threads 10 --repeatMaskerAnnotation repeat_masker_annotation_file
    }
    runtime {
     docker: docker_image
    }
    output {
     File pmd_all = "pmd-all.csv"
     File umr_lmr_all = "umr-lmr-all.csv"
     File pmd_segments = "pmd-segments.csv"
     File umr_lmr_with_pmd = "./LMRUMRwithPMD/umr-lmr.csv"
     File umr_lmr_without_pmd = "./LMRUMRwithoutPMD/umr-lmr.csv"
    }
}
