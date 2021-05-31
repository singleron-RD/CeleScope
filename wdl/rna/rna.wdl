version 1.0

import "../tools/common.wdl" as step_common
import "star.wdl" as step_star
import "../tools/featureCounts.wdl" as step_featureCounts
import "count.wdl" as step_count
import "analysis.wdl" as step_analysis


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

        String docker_use

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
        Int? cpu_analysis
        Int? mem_analysis
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
            docker_use = docker_use,
    }


    call step_star.star {
            
        input:
            sample_name = sample_name,
            cutadapt_out_fq = run_common.cutadapt_out_fq,
            genomeDir = genomeDir,
            in_data = run_common.out_data,
            cpu_star = cpu_star,
            mem_star = mem_star,
            docker_use = docker_use
    } 

    call step_featureCounts.featureCounts {            
        input:
            sample_name = sample_name,
            in_bam = star.out_bam,
            gtf_type = gtf_type,
            genomeDir = genomeDir,
            in_data = star.out_data,
            cpu_featureCounts = cpu_featureCounts,
            mem_featureCounts = mem_featureCounts,
            docker_use = docker_use,
    } 

    call step_count.count {
        input:
            sample_name = sample_name,
            in_bam = featureCounts.out_bam,
            genomeDir = genomeDir,
            in_data = featureCounts.out_data,
            mem_on_bam = featureCounts.mem_on_bam,
            cpu_count = cpu_count,
            mem_count = mem_count,
            docker_use = docker_use,
    }

    call step_analysis.analysis {
        input:
            sample_name = sample_name,
            in_matrix = count.out_matrix,
            in_data = count.out_data,
            mem_on_mtx = count.mem_on_mtx,
            genomeDir = genomeDir,

            cpu_analysis = cpu_analysis,
            mem_analysis = mem_analysis,
            docker_use = docker_use,
    }

    output {
        File bam = featureCounts.out_bam
        File h5ad = analysis.out_h5ad
        File report = analysis.out_report
    }


}