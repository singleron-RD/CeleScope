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
from celescope.tools.utils import merge_report, generate_sjm


def parse_map(mapfile):
    fq_dict = defaultdict(list)
    match_dict = defaultdict(list)
    with open(mapfile) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                continue
            tmp = line.split()
            library_id = tmp[0]
            library_path = tmp[1]
            sample_name = tmp[2]
            match_dir = tmp[3]

            try:
                pattern1_1 = library_path + '/' + library_id + '*' + '_1.fq.gz'
                pattern1_2 = library_path + '/' + library_id + '*' + 'R1_*.fastq.gz'
                pattern2_1 = library_path + '/' + library_id + '*' + '_2.fq.gz'
                pattern2_2 = library_path + '/' + library_id + '*' + 'R2_*.fastq.gz'
                fq1 = (glob.glob(pattern1_1) + glob.glob(pattern1_2))[0]
                fq2 = (glob.glob(pattern2_1) + glob.glob(pattern2_2))[0]
            except IndexError as e:
                sys.exit("Mapfile Error:" + str(e))

            assert os.path.exists(fq1), '%s not exists!' % (fq1)
            assert os.path.exists(fq2), '%s not exists!' % (fq2)
            if sample_name in fq_dict:
                fq_dict[sample_name][0].append(fq1)
                fq_dict[sample_name][1].append(fq2)
            else:
                fq_dict[sample_name] = [[fq1], [fq2]]
            match_dict[sample_name] = match_dir

    for sample_name in fq_dict:
        fq_dict[sample_name][0] = ",".join(fq_dict[sample_name][0])
        fq_dict[sample_name][1] = ",".join(fq_dict[sample_name][1])

    return fq_dict, match_dict


def main():

    parser = argparse.ArgumentParser('CeleScope virus multi-sample')
    #parser.add_argument('--mod', help='mod, sjm or shell', choices=['sjm', 'shell'], default='sjm')
    parser.add_argument(
        '--mapfile',
        help='mapfile, 3 columns, "LibName\\tDataDir\\tSampleName"',
        required=True)
    parser.add_argument('--chemistry', choices=['scopeV2.0.0', 'scopeV2.0.1',
                                                'scopeV2.1.0', 'scopeV2.1.1'], help='chemistry version')
    parser.add_argument('--whitelist', help='cellbarcode list')
    parser.add_argument('--linker', help='linker')
    parser.add_argument('--pattern', help='read1 pattern')
    parser.add_argument('--outdir', help='output dir', default="./")
    parser.add_argument(
        '--adapt',
        action='append',
        help='adapter sequence',
        default=[
            'polyT=A{15}',
            'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
    parser.add_argument(
        '--minimum-length',
        dest='minimum_length',
        help='minimum_length',
        default=20)
    parser.add_argument(
        '--nextseq-trim',
        dest='nextseq_trim',
        help='nextseq_trim',
        default=20)
    parser.add_argument(
        '--overlap',
        help='minimum overlap length, default=5',
        default=5)
    parser.add_argument(
        '--lowQual',
        type=int,
        help='max phred of base as lowQual',
        default=0)
    parser.add_argument(
        '--lowNum',
        type=int,
        help='max number with lowQual allowed',
        default=2)
    parser.add_argument('--starMem', help='starMem', default=30)
    parser.add_argument('--genomeDir', help='genome index dir', required=True)
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument(
        '--virus_genomeDir',
        help='virus_genomeDir',
        required=True)
    args = vars(parser.parse_args())

    fq_dict, match_dict = parse_map(args['mapfile'])

    # 链接数据
    raw_dir = args['outdir'] + '/data_give/rawdata'
    os.system('mkdir -p %s' % (raw_dir))
    with open(raw_dir + '/ln.sh', 'w') as fh:
        fh.write('cd %s\n' % (raw_dir))
        for s, arr in fq_dict.items():
            fh.write('ln -sf %s %s\n' % (arr[0], s + '_1.fq.gz'))
            fh.write('ln -sf %s %s\n' % (arr[1], s + '_2.fq.gz'))
    #os.system('sh %s'%(raw_dir+'/ln.sh'))

    logdir = args['outdir'] + '/log'
    os.system('mkdir -p %s' % (logdir))
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    conda = __CONDA__
    app = 'celescope'
    thread = args['thread']
    chemistry = args['chemistry']
    genomeDir = args['genomeDir']
    virus_genomeDir = args['virus_genomeDir']
    pattern = args['pattern']
    whitelist = args['whitelist']
    linker = args['linker']
    lowQual = args['lowQual']
    lowNum = args['lowNum']
    starMem = args['starMem']
    assay = "capture_virus"

    basedir = args['outdir']
    steps = [
        'sample',
        'barcode',
        'cutadapt',
        "STAR_virus",
        "count_capture_virus",
        "analysis_capture_virus"]

    for sample in fq_dict:
        outdir_dic = {}
        index = 0
        for step in steps:
            outdir = f"{basedir}/{sample}/{index:02d}.{step}"
            outdir_dic.update({step: outdir})
            index += 1

        # sample
        step = "sample"
        cmd = f'''source activate {conda};  {app} {assay} {step} --chemistry {chemistry}
        --sample {sample} --outdir {outdir_dic[step]} --genomeDir {genomeDir} --assay {assay}'''
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}')
        last_step = step

        # barcode
        arr = fq_dict[sample]
        step = "barcode"
        cmd = f'''source activate {conda}; {app} {assay} {step}  --fq1 {arr[0]} --fq2 {arr[1]} --chemistry {chemistry}
            --pattern {pattern} --whitelist {whitelist} --linker {linker} --sample {sample} --lowQual {lowQual}
            --lowNum {lowNum} --outdir {outdir_dic[step]} --thread {thread} --assay {assay}'''
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', m=5, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        last_step = step

        # adapt
        step = "cutadapt"
        fq = f'{outdir_dic["barcode"]}/{sample}_2.fq.gz'
        cmd = f'''source activate {conda};  {app} {assay} {step}  --fq {fq} --sample {sample} --outdir
            {outdir_dic[step]} --assay {assay}'''
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        last_step = step

        # STAR_virus
        step = 'STAR_virus'
        input_read = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = f'''source activate {conda};  {app} {assay} {step}  --input_read {input_read} --sample {sample}
        --virus_genomeDir {virus_genomeDir} --thread {thread} --outdir {outdir_dic[step]} --assay {assay}'''
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', m=starMem, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        last_step = step

        # count_capture_virus
        step = 'count_capture_virus'
        virus_bam = f'{outdir_dic["STAR_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
        cmd = f'''source activate {conda};  {app} {assay} {step}  --virus_bam {virus_bam} --sample {sample}
            --outdir {outdir_dic[step]} --match_dir {match_dict[sample]} --assay {assay}'''
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', m=20, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        last_step = step

        # analysis_capture_virus
        step = 'analysis_capture_virus'
        virus_file = f'{outdir_dic["count_capture_virus"]}/{sample}_virus_UMI_count.tsv'
        cmd = f'''source activate {conda};  {app} {assay} {step}  --match_dir {match_dict[sample]} --sample {sample}
            --outdir {outdir_dic[step]} --virus_file {virus_file} --assay {assay}'''
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', m=15, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        last_step = step

    # merged report
    merge_report(fq_dict, steps, last_step, sjm_cmd, sjm_order, logdir, conda)


if __name__ == '__main__':
    main()
