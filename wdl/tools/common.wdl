version 1.0

workflow run_common {
    input {
        String sample_name
        Array[File] raw_fq1s
        Array[File] raw_fq2s
        String? pattern = "None"
        String? whitelist = "None"
        String? linker = "None"
        Int lowqual = 0
        Int lownum = 2
        Int overlap = 10
        Int min_length = 20
        Int insert = 150

        Int? cpu_sample
        Int? mem_sample
        Int? cpu_barcode
        Int? mem_barcode
        Int? cpu_cutadapt
        Int? mem_cutadapt
    }

    Int cpu_default = 1
    Int mem_default = 1

    Int runtime_cpu_sample = select_first([cpu_sample, cpu_default])
    Int runtime_mem_sample = select_first([mem_sample, mem_default])

    Int runtime_cpu_barcode = select_first([cpu_barcode, cpu_default])
    Int runtime_mem_barcode = select_first([mem_barcode, mem_default])

    Int runtime_cpu_cutadapt = select_first([cpu_cutadapt, cpu_default])
    Int runtime_mem_cutadapt = select_first([mem_cutadapt, mem_default])

    call sample {
        input:
            sample_name = sample_name,
            raw_fq1s = raw_fq1s,
            runtime_cpu_sample = runtime_cpu_sample,
            runtime_mem_sample = runtime_mem_sample,
    }


    call barcode {
        input:
        sample_name = sample_name,
        raw_fq1s = raw_fq1s,
        raw_fq2s = raw_fq2s,
        pattern = pattern,
        whitelist = whitelist,
        linker = linker,
        lowqual = lowqual,
        lownum = lownum,
        in_data = sample.out_data,
        runtime_cpu_barcode = runtime_cpu_barcode,
        runtime_mem_barcode = runtime_mem_barcode,
    }

    call cutadapt {
        input:
        sample_name = sample_name,
        valid_fq = barcode.valid_fq,
        overlap = overlap,
        min_length = min_length,
        insert = insert,
        in_data = barcode.out_data,
        runtime_cpu_cutadapt = runtime_cpu_cutadapt,
        runtime_mem_cutadapt = runtime_mem_cutadapt,
    }

    output {
        File out_data = cutadapt.out_data
        File cutadapt_out_fq = cutadapt.cutadapt_out_fq     
    }
}

task sample {
    input {
        String sample_name
        Array[File] raw_fq1s
        Int runtime_cpu_sample
        Int runtime_mem_sample
    }

    runtime {
        cpu: runtime_cpu_sample
        memory: runtime_mem_sample + "GiB"
    }

    command {
        set -euo pipefail
        celescope rna sample \
        --outdir "00.sample" \
        --sample "~{sample_name}" \
        --assay "rna" \
        --chemistry "auto" \
        --fq1 "~{sep="," raw_fq1s}"
        ls -alh
    }

    output {
        File out_data = ".data.json"
    }
}

task barcode {
    input {
        String sample_name
        Array[File] raw_fq1s
        Array[File] raw_fq2s
        String? pattern
        String? whitelist
        String? linker
        Int lowqual
        Int lownum
        File in_data
        Int runtime_cpu_barcode
        Int runtime_mem_barcode
    }

    runtime {
        cpu: runtime_cpu_barcode
        memory: runtime_mem_barcode + "GiB"
    }


    command {
        set -euo pipefail
        mv "~{in_data}" ".data.json"
        celescope rna barcode \
        --outdir "01.barcode" \
        --sample "~{sample_name}" \
        --assay "rna" \
        --chemistry "auto" \
        --fq1 "~{sep="," raw_fq1s}" \
        --fq2 "~{sep="," raw_fq2s}" \
        --pattern "~{default="None" pattern}" \
        --whitelist "~{default="None" whitelist}" \
        --linker "~{default="None" linker}" \
        --lowQual "~{lowqual}" \
        --lowNum "~{lownum}"
    }


    output {
        File out_data = ".data.json"
        File valid_fq = "01.barcode/~{sample_name}_2.fq"
    }
}

task cutadapt {

    input {
        String sample_name
        File valid_fq
        Int overlap
        Int min_length
        Int insert
        File in_data
        Int runtime_cpu_cutadapt
        Int runtime_mem_cutadapt
    }

    runtime {
        cpu: runtime_cpu_cutadapt
        memory: runtime_mem_cutadapt + "GiB"
    }

    command {
        set -euo pipefail
        mv "~{in_data}" ".data.json"
        celescope rna cutadapt \
        --outdir "02.cutadapt" \
        --sample "~{sample_name}" \
        --assay "rna" \
        --fq "~{valid_fq}" \
        --overlap "~{overlap}" \
        --minimum_length "~{min_length}" \
        --insert "~{insert}" \
        --thread "~{runtime_cpu_cutadapt}"
    }

    output {
        File out_data = ".data.json"
        File cutadapt_out_fq = "02.cutadapt/~{sample_name}_clean_2.fq"
    }
}