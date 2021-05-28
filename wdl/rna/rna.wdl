version 1.0

import "../tools/common.wdl" as step_common

workflow run_rna {

    input {
        String sample_name
        Array[File] raw_fq1s
        Array[File] raw_fq2s
    }

    call step_common.run_common {
        input:
            sample_name = sample_name,
            raw_fq1s = raw_fq1s,
            raw_fq2s = raw_fq2s
    }

    output {
        File data_json = run_common.data_json
        File cutadapt_out_fq = run_common.cutadapt_out_fq
    }

}