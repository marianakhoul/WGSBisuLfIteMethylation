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
     File bedgraph_ratio = "${methylation_dir}/${sample_name}_CpG_ratio.bedGraph"
    }
}

task metilene_input {
    
    String metilene_dir
    File bedgraph_ratio
    String methylation_dir
    String docker_image
    
    command {
     bedtools unionbedg -filler NA -i ${methylation_dir}/${bedgraph_ratio} > ${metilene_dir}/metilene-input.bedGraph
    }
    runtime {
     docker: docker_image
    }
    output {
     File metilene_input_file = "${metilene_dir}/metilene-input.bedGraph"
    }
}


task metilene {

    String metilene_dir
    File log
    String min_cpg
    String min_diff
    String threads
    File metilene_input_file
    String docker_image
    
    command {
     metilene -m ${min_cpg} -d ${min_diff} -t ${threads} ${metilene_dir}/${metilene_input_file} > ${metilene_dir}/dmrs.csv 2> ${log}
    }
    runtime {
     docker: docker_image
    }
    output {
     File metilene_output = "${metilene_dir}/dmrs.csv"
    }
}


task camel_index {

    String docker_image
    String camel_modules_path
    File ref_fasta
    String reference_fasta
    
    command {
     python ${camel_modules_path}/index.py ${ref_fasta} ${reference_fasta}.h5
    }
    runtime {
     docker: docker_image
    }
    output {
     File reference_output = "${reference_fasta}.h5"
    }
}

task camel_call {

    String camel_modules_path
    File reference_output
    File log
    String camel_dir
    String sample_name
    File input_bam
    File bam_bai
    String docker_image
    
    command {
     python ${camel_modules_path}/call.py ${input_bam} ${reference_output} ${camel_dir}/${sample_name}.h5 2> ${log}
    }
    runtime {
     docker: docker_image
    }
    output {
     File camel_call_output = "${camel_dir}/{sample_name}.h5"
    }
}

task camel_dmr {

    File reference_output
    String camel_dir
    String docker_image
    File log
    String min_diff
    String min_cpg
    String min_cov
    
    command {
     python ${camel_modules_path}/dmr.py ${reference_output} \
     --case {input.case} \
     --control {input.control} \
     --min_diff ${min_diff} \
     --min_cpg ${min_cpg} \
     --min_cov ${min_cov} > ${camel_dir}/dmrs.csv 2> ${log}
    }
    runtime {
     docker: docker_image
    }
    output {
     File camel_call_output = "${camel_dir}/dmrs.csv"
    }

}
