version 1.0

import "structs.wdl"

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

        Runtime runtime_sample
        Runtime runtime_barcode
        Runtime runtime_cutadapt
    }

    call sample {
        input:
            sample_name = sample_name,
            raw_fq1s = raw_fq1s,

            runtime_sample = runtime_sample,
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
        runtime_barcode = runtime_barcode,
    }

    call cutadapt {
        input:
        sample_name = sample_name,
        valid_fq = barcode.valid_fq,
        overlap = overlap,
        min_length = min_length,
        insert = insert,
        in_data = barcode.out_data,
        runtime_cutadapt = runtime_cutadapt,
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
        Runtime runtime_sample

    }

    runtime {
        cpu: runtime_sample.cpu
        memory: runtime_sample.memory_gb + "GiB"
        docker: runtime_sample.docker
        queue: runtime_sample.queue
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
        Runtime runtime_barcode
    }

    runtime {
        cpu: runtime_barcode.cpu
        memory: runtime_barcode.memory_gb + "GiB"
        docker: runtime_barcode.docker
        queue: runtime_barcode.queue
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
        Runtime runtime_cutadapt
    }

    runtime {
        cpu: runtime_cutadapt.cpu
        memory: runtime_cutadapt.memory_gb + "GiB"
        docker: runtime_cutadapt.docker
        queue: runtime_cutadapt.queue
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
        --thread "~{runtime_cutadapt.cpu}"
    }

    output {
        File out_data = ".data.json"
        File cutadapt_out_fq = "02.cutadapt/~{sample_name}_clean_2.fq"
    }
}