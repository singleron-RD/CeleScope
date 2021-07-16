version 1.0

import "structs.wdl"

task featureCounts {
    input {
        String sample_name
        File in_bam
        String gtf_type
        String genomeDir
        File in_data
        Runtime runtime_featureCounts
    }

    runtime {
        cpu: runtime_featureCounts.cpu
        memory: runtime_featureCounts.memory_gb + "GiB"
        docker: runtime_featureCounts.docker
        queue: runtime_featureCounts.queue
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
        --thread ~{runtime_featureCounts.cpu}
    }

    output {
        File out_data = ".data.json"
        File out_bam = "04.featureCounts/~{sample_name}_name_sorted.bam"
        Int mem_on_bam = ceil(size(out_bam, "G") * 3 + 3)
    }
}


