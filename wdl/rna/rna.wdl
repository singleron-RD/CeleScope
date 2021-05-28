version 1.0

import "../tools/common.wdl" as step_common
import "star.wdl" as step_star

workflow run_rna {

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
        String genomeDir
        String gtf_type = "exon"

        Int? cpu_sample
        Int? mem_sample
        Int? cpu_barcode
        Int? mem_barcode
        Int? cpu_cutadapt
        Int? mem_cutadapt
        Int? cpu_star
        Int? mem_star
        Int? cpu_featureCounts
        Int? mem_featureCounts
        Int? cpu_count
        Int? mem_count
    }


    call step_common.run_common {
        input:
            sample_name = sample_name,
            raw_fq1s = raw_fq1s,
            raw_fq2s = raw_fq2s,
            cpu_sample = cpu_sample,
            mem_sample = mem_sample,
            cpu_barcode = cpu_barcode,
            mem_barcode = mem_barcode,
            cpu_cutadapt = cpu_cutadapt,
            mem_cutadapt = mem_cutadapt,      
    }

    Int runtime_cpu_star = select_first([cpu_star, 6])
    Int runtime_mem_star = select_first([mem_star, 30])

    call step_star.star {
            
        input:
            sample_name = sample_name,
            cutadapt_out_fq = run_common.cutadapt_out_fq,
            genomeDir = genomeDir,
            in_data = run_common.data_json,
            runtime_cpu_star = runtime_cpu_star,
            runtime_mem_star = runtime_mem_star,
    } 

    output {
        File data_json = star.data_json
    }

}