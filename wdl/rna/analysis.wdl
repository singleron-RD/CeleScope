version 1.0

task analysis {
    input {
        String sample_name
        File in_matrix
        File in_data
        Int mem_on_mtx
        String genomeDir

        Int? cpu_analysis
        Int? mem_analysis
    }

    Int runtime_cpu_analysis = select_first([cpu_analysis, 1])
    Int runtime_mem_analysis = select_first([mem_analysis, 8])

    runtime {
        cpu: runtime_cpu_analysis
        memory:  if mem_on_mtx > runtime_mem_analysis then mem_on_mtx + "GiB" else runtime_mem_analysis + "GiB"
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
        --thread "~{runtime_cpu_analysis}"
    }

    output {
        File out_data = ".data.json"
        File out_h5ad = "06.analysis/~{sample_name}.h5"
        File out_report = "~{sample_name}_report.html"
    }
}
