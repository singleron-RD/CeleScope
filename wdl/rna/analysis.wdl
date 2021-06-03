version 1.0

import "../tools/structs.wdl"

task analysis {
    input {
        String sample_name
        File in_matrix
        File in_data
        Int mem_on_mtx
        String genomeDir
        Runtime runtime_analysis
    }


    runtime {
        cpu: runtime_analysis.cpu
        memory:  if mem_on_mtx > runtime_analysis.memory_gb then mem_on_mtx + "GiB" else runtime_analysis.memory_gb + "GiB"
        docker: runtime_analysis.docker
        queue: runtime_analysis.queue
    }

    command {
        set -euo pipefail
        mkdir -p "05.count/~{sample_name}_matrix_10X"
        tar xf "~{in_matrix}" --directory="05.count/~{sample_name}_matrix_10X"
        mv "~{in_data}" ".data.json"
        celescope rna analysis \
        --outdir "06.analysis" --sample "~{sample_name}" --assay rna \
        --matrix_file "05.count/~{sample_name}_matrix_10X" \
        --genomeDir "~{genomeDir}" \
        --thread "~{runtime_analysis.cpu}"
    }

    output {
        File out_data = ".data.json"
        File out_h5ad = "06.analysis/~{sample_name}.h5"
        File out_report = "~{sample_name}_report.html"
    }
}
