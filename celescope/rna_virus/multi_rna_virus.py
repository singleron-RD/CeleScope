#!/bin/env python
# coding=utf8

import os
import glob
import sys
import argparse
import re
import logging
from collections import defaultdict
from celescope.__init__ import __CONDA__
from celescope.rna_virus.__init__ import __STEPS__, __ASSAY__
from celescope.tools.utils import merge_report, generate_sjm
from celescope.tools.utils import parse_map_col4, multi_opts, link_data


def main():

    # init
    assay = __ASSAY__
    steps = __STEPS__
    conda = __CONDA__
    app = 'celescope'

    # parser
    parser = multi_opts(assay)
    parser.add_argument('--starMem', help='starMem', default=30)
    parser.add_argument('--genomeDir', help='genome index dir', required=True)
    parser.add_argument(
        '--gtf_type',
        help='Specify attribute type in GTF annotation, default=exon',
        default='exon')
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument(
        '--virus_genomeDir',
        help='virus_genomeDir',
        required=True)
    args = parser.parse_args()

    fq_dict, cells_dict = parse_map_col4(args.mapfile, 'auto')

    # read args
    outdir = args.outdir
    chemistry = args.chemistry
    pattern = args.pattern
    whitelist = args.whitelist
    linker = args.linker
    lowQual = args.lowQual
    lowNum = args.lowNum
    mod = args.mod
    rm_files = args.rm_files

    # parse mapfile
    fq_dict, cells_dict = parse_map_col4(args.mapfile, "auto")

    # link
    link_data(outdir, fq_dict)

    # custom args
    thread = args.thread
    genomeDir = args.genomeDir
    starMem = args.starMem
    gtf_type = args.gtf_type
    virus_genomeDir = args.virus_genomeDir

    # mk log dir
    logdir = outdir + '/log'
    os.system('mkdir -p %s' % (logdir))

    # script init
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    shell_dict = defaultdict(str)

    # outdir dict
    for sample in fq_dict:
        outdir_dic = {}
        index = 0
        for step in steps:
            step_outdir = f"{outdir}/{sample}/{index:02d}.{step}"
            outdir_dic.update({step: step_outdir})
            index += 1

        # sample
        step = "sample"
        cmd = (
            f'{app} {assay} {step} '
            f'--chemistry {chemistry} '
            f'--sample {sample} --outdir {outdir_dic[step]} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda)
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # barcode
        arr = fq_dict[sample]
        step = "barcode"
        cmd = (
            f'{app} {assay} {step} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} --chemistry {chemistry} '
            f'--pattern {pattern} --whitelist {whitelist} --linker {linker} '
            f'--sample {sample} --lowQual {lowQual} --thread {thread} '
            f'--lowNum {lowNum} --outdir {outdir_dic[step]} --assay {assay} '

        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # adapt
        step = "cutadapt"
        fq = f'{outdir_dic["barcode"]}/{sample}_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--fq {fq} --sample {sample} --outdir '
            f'{outdir_dic[step]} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # STAR
        step = 'STAR'
        fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--fq {fq} --sample {sample} '
            f'--outdir {outdir_dic[step]} --assay {assay} '
            f'--genomeDir {genomeDir} --thread {thread} '
            f'--out_unmapped '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=starMem, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # STAR_virus
        step = 'STAR_virus'
        input_read = f'{outdir_dic["STAR"]}/{sample}_Unmapped.out.mate1'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} --outdir {outdir_dic[step]} --assay {assay} '
            f'--input_read {input_read} '
            f'--virus_genomeDir {virus_genomeDir} '
            f'--thread {thread} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=starMem, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # featureCounts
        step = 'featureCounts'
        input = f'{outdir_dic["STAR"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{app} {assay} {step} '
            f'--input {input} --gtf_type {gtf_type} '
            f'--sample {sample} --thread {thread} --outdir {outdir_dic[step]} '
            f'--genomeDir {genomeDir} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=8, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # count
        step = 'count'
        bam = f'{outdir_dic["featureCounts"]}/{sample}_name_sorted.bam'
        cmd = (
            f'{app} {assay} {step} '
            f'--bam {bam} --sample {sample} --cells {cells_dict[sample]} '
            f'--outdir {outdir_dic[step]} --assay {assay} '
            f'--genomeDir {genomeDir}'
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=8, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # count_virus
        step = 'count_virus'
        virus_bam = f'{outdir_dic["STAR_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
        barcode_file = f'{outdir_dic["count"]}/{sample}_matrix_10X/barcodes.tsv'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} --outdir {outdir_dic[step]} --assay {assay} '
            f'--virus_bam {virus_bam} '
            f'--barcode_file {barcode_file} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=20, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

        # analysis_rna_virus
        step = 'analysis_rna_virus'
        matrix_file = f'{outdir_dic["count"]}/{sample}_matrix.tsv.gz'
        virus_file = f'{outdir_dic["count_virus"]}/{sample}_virus_UMI_count.tsv'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} --outdir {outdir_dic[step]} --assay {assay} '
            f'--matrix_file {matrix_file} '
            f'--virus_file {virus_file} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=15, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell_dict[sample] += cmd + '\n'
        last_step = step

    # merged report
    if mod == 'sjm':
        step = 'merge_report'
        merge_report(
            fq_dict,
            steps,
            last_step,
            sjm_cmd,
            sjm_order,
            logdir,
            conda,
            outdir,
            rm_files,
        )
    if mod == 'shell':
        os.system('mkdir -p ./shell/')
        for sample in shell_dict:
            with open(f'./shell/{sample}.sh', 'w') as f:
                f.write(shell_dict[sample])


if __name__ == '__main__':
    main()
