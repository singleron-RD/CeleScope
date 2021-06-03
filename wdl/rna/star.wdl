version 1.0

import "../tools/structs.wdl"

task star {
    input {
        String sample_name
        File cutadapt_out_fq
        String genomeDir
        File in_data
        Runtime runtime_star

    }

    runtime {
        cpu: runtime_star.cpu
        memory: runtime_star.memory_gb + "GiB"
        docker: runtime_star.docker
        queue: runtime_star.queue
    }

    command {
        set -euo pipefail
        mv "~{in_data}" ".data.json"
        celescope rna star \
        --outdir "03.STAR" \
        --sample "~{sample_name}" \
        --assay rna --fq "~{cutadapt_out_fq}" \
        --genomeDir "~{genomeDir}" \
        --thread ~{runtime_star.cpu} \
        --starMem ~{runtime_star.memory_gb} \
        --outFilterMatchNmin 0
    }

    output {
        File out_data = ".data.json"
        File out_bam = "03.STAR/~{sample_name}_Aligned.sortedByCoord.out.bam"
    }
}