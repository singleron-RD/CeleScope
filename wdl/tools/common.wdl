version 1.0

struct RuntimeAttr {
    Int? cpu
    Int? memory_gb
}

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

        RuntimeAttr? runtime_attr_sample
        RuntimeAttr? runtime_attr_barcode
        RuntimeAttr? runtime_attr_cutadapt
    }

    call sample {
        input:
            sample_name = sample_name,
            raw_fq1s = raw_fq1s
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
        in_data = sample.data_json,

        runtime_attr_override = runtime_attr_barcode
    }

    call cutadapt {
        input:
        sample_name = sample_name,
        valid_fq = barcode.valid_fq,
        overlap = overlap,
        min_length = min_length,
        insert = insert,
        in_data = barcode.data_json,

        runtime_attr_override = runtime_attr_cutadapt
    }

    output {
        File data_json = cutadapt.data_json
        File cutadapt_out_fq = cutadapt.cutadapt_out_fq     
    }
}

task sample {
    input {
        String sample_name
        Array[File] raw_fq1s
        RuntimeAttr? runtime_attr_override
    }
    command {
        set -euo pipefail
        celescope rna sample \
        --outdir "00.sample" --sample "~{sample_name}" --assay "rna" \
        --chemistry "auto" --fq1 "~{sep="," raw_fq1s}"
        ls -alh
    }
    RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1,
}
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])

    runtime {
        cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
        memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    }
    output {
        File data_json = ".data.json"
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

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_attr_default = object {
        cpu: 1,
        memory_gb: 1,
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])

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
        --thread ~{runtime_attr.cpu} \
        --lowNum "~{lownum}"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
        memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    }
    output {
        File data_json = ".data.json"
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

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_attr_default = object {
        cpu: 1,
        memory_gb: 1,
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])

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
        --thread "~{runtime_attr.cpu}"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
        memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    }
    output {
        File data_json = ".data.json"
        File cutadapt_out_fq = "02.cutadapt/~{sample_name}_clean_2.fq"
    }
}