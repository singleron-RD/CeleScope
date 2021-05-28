version 1.0

task featureCounts {
    input {
        String sample_name
        File in_bam
        String gtf_type
        String genomeDir
        File in_data
        Int? cpu_featureCounts
        Int? mem_featureCounts
    }

    Int runtime_cpu_featureCounts = select_first([cpu_featureCounts, 6])
    Int runtime_mem_featureCounts = select_first([mem_featureCounts, 2])

    runtime {
        cpu: runtime_cpu_featureCounts
        memory: runtime_mem_featureCounts + "GiB"
    }

    command {
        set -euo pipefail
        mv "~{in_data}" ".data.json"
        celescope rna featureCounts \
        --outdir "04.featureCounts" \
        --sample "~{sample_name}" \
        --assay "rna" \
        --input "~{in_bam}" \
        --gtf_type "~{gtf_type}" \
        --genomeDir "~{genomeDir}" \
        --thread ~{runtime_cpu_featureCounts}
    }

    output {
        File out_data = ".data.json"
        File out_bam = "04.featureCounts/~{sample_name}_name_sorted.bam"
        Int mem_on_bam = ceil(size(out_bam, "G") * 3 + 3)
    }
}


