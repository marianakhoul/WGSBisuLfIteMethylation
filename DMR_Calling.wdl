version 1.1

#
# Task associated with DMR calling for the whole Genome BisuLfIte 
# sequencing Methylation analysis Pipeline.
#

task methyl_dackel {
    input {    
        String docker_image
        String sample_name
        File ref_fasta_index
        File ref_fasta
        File bam_index
        File bam_file
    }
    
    command {
        MethylDackel extract --mergeContext \
            -o ${sample_name} ${ref_fasta} ${bam_file}
    }

    runtime {
        docker: docker_image
    }

    output {
        File methyl_dackel_output = "${sample_name}_CpG.bedGraph"
    }
}

task bedgraph_to_methylation_ratio {
    input {
        File methyl_dackel_output
        String sample_name
        String docker_image
    }

    command {
        Rscript /usr/local/bin/transformBedGraph.R
    }

    runtime {
        docker: docker_image
    }

    output {
        File bedgraph_ratio = "${sample_name}_CpG_ratio.bedGraph"
    }
}

task metilene_input {
    input { 
        File bedgraph_ratio
        String docker_image
    } 
    
    command {
        bedtools unionbedg -filler NA \
            -i ${bedgraph_ratio} > metilene-input.bedGraph
    }

    runtime {
        docker: docker_image
    }

    output {
        File metilene_input_file = "metilene-input.bedGraph"
    }
}

task metilene {
    input {
        Int min_cpg
        Float min_diff
        Int threads
        File metilene_input_file
        String docker_image
    }
    
    command {
        metilene -m ${min_cpg} -d ${min_diff} -t ${threads} \
            ${metilene_input_file} > dmrs.csv
    }

    runtime {
        docker: docker_image
    }

    output {
        File metilene_output = "dmrs.csv"
    }
}

task bsseq {
    input {    
        String docker_image
        File bedgraph_ratio
        String wg_blimp_R_script_path
        Int io_threads
    }
    
    command {
        R ${wg_blimp_R_script_path}/bsseq.R
    }

    runtime {
        docker: docker_image
    }

    output {
        File rdata_file = "bsseq.Rdata"
        File csv_file = "dmrs.csv"
        File pdf_file = "top100.pdf"
    }
}
