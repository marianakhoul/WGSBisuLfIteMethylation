## Task associated with DMR calling for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task methyl_dackel {
    
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_fasta_index
    File ref_fasta
    File input_bam
    String sample_name
    String methylation_dir
    
    
    command {
     MethylDackel extract --mergeContext -o ${methylation_dir}/${sample_name} ${ref_fasta} ${input_bam}
    }
    runtime {
     docker: docker_image
    }
    output {
     File methyl_dackel_output = "${methylation_dir}/${sample_name}_CpG.bedGraph"
    }

}
