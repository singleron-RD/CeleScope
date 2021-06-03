version 1.0

import "tools/common.wdl" as step_common
import "rna/star.wdl" as step_star
import "tools/featureCounts.wdl" as step_featureCounts
import "rna/count.wdl" as step_count
import "rna/analysis.wdl" as step_analysis
import "tools/structs.wdl"


workflow rna {

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

        Runtime runtime_sample
        Runtime runtime_barcode
        Runtime runtime_cutadapt
        Runtime runtime_star
        Runtime runtime_featureCounts
        Runtime runtime_count
        Runtime runtime_analysis
    }


    call step_common.run_common {
        input:
            sample_name = sample_name,
            raw_fq1s = raw_fq1s,
            raw_fq2s = raw_fq2s,

            runtime_sample = runtime_sample,
            runtime_barcode = runtime_barcode,
            runtime_cutadapt = runtime_cutadapt,
    }


    call step_star.star {
            
        input:
            sample_name = sample_name,
            cutadapt_out_fq = run_common.cutadapt_out_fq,
            genomeDir = genomeDir,
            in_data = run_common.out_data,

            runtime_star = runtime_star,
    } 

    call step_featureCounts.featureCounts {            
        input:
            sample_name = sample_name,
            in_bam = star.out_bam,
            gtf_type = gtf_type,
            genomeDir = genomeDir,
            in_data = star.out_data,

            runtime_featureCounts = runtime_featureCounts,
    } 

    call step_count.count {
        input:
            sample_name = sample_name,
            in_bam = featureCounts.out_bam,
            genomeDir = genomeDir,
            in_data = featureCounts.out_data,
            mem_on_bam = featureCounts.mem_on_bam,
            runtime_count = runtime_count,
    }

    call step_analysis.analysis {
        input:
            sample_name = sample_name,
            in_matrix = count.out_matrix,
            in_data = count.out_data,
            mem_on_mtx = count.mem_on_mtx,
            genomeDir = genomeDir,

            runtime_analysis = runtime_analysis
    }

    output {
        File bam = featureCounts.out_bam
        File h5ad = analysis.out_h5ad
        File report = analysis.out_report
    }


}