## Task associated with DMR calling for the whole Genome BisuLfIte sequencing Methylation analysis Pipeline.
##
##
## Task Description

task methyl_dackel {
    
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
    String methylation_dir
    File log
    
    command {
     MethylDackel extract --mergeContext -o ${methylation_dir}/${sample_name} ${ref_fasta} ${bam_file}
    }
    runtime {
     docker: docker_image
    }
    output {
     File methyl_dackel_output = "${methylation_dir}/${sample_name}_CpG.bedGraph"
    }

}

task bedgraph_to_methylation_ratio {
    
    String methylation_dir
    File log
    File methyl_dackel_output
    String sample_name
    String wg_blimp_R_script_path
    String docker_image
    
    command {
     R ${wg_blimp_R_script_path}/transformBedGraph.R
    }
    runtime {
     docker: docker_image
    }
    output {
     File bedgraph_to_methylation_ratio = "${methylation_dir}/${sample_name}_CpG_ratio.bedGraph"
    }
}

task metilene_input {
    
    String metilene_dir
    File bedgraph_to_methylation_ratio
    String methylation_dir
    String docker_image
    
    command {
     bedtools unionbedg -filler NA -header -names {g1_header} {g2_header} -i ${methylation_dir}/${bedgraph_to_methylation_ratio} > ${metilene_dir}/${metilene_input}
    }
    runtime {
     docker: docker_image
    }
    output {
     File metilene_input = "${metilene_dir}/metilene-input.bedGraph"
    }
}
