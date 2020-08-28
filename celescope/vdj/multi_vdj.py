import os
import glob
import sys
import argparse
import re
from collections import defaultdict
from celescope.__init__ import __CONDA__
from celescope.vdj.__init__ import __STEPS__, __ASSAY__
from celescope.tools.utils import merge_report, generate_sjm, parse_map_col4


def main():

    parser = argparse.ArgumentParser('CeleScope vdj multi-sample')
    parser.add_argument(
        '--mapfile', help='mapfile, 3 columns, "LibName\\tDataDir\\tSampleName"', required=True)
    parser.add_argument('--mod', help='mod, sjm or shell', choices=['sjm', 'shell'], default='sjm')
    parser.add_argument('--chemistry', choices=['scopeV2.0.0', 'scopeV2.0.1',
                                                'scopeV2.1.0', 'scopeV2.1.1'], help='chemistry version')
    parser.add_argument('--whitelist', help='cellbarcode list')
    parser.add_argument('--linker', help='linker')
    parser.add_argument('--pattern', help='read1 pattern')
    parser.add_argument('--outdir', help='output dir', default="./")
    parser.add_argument('--adapt', action='append', help='adapter sequence',
                        default=['polyT=A{15}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
    parser.add_argument('--minimum-length', dest='minimum_length',
                        help='minimum_length', default=20)
    parser.add_argument('--nextseq-trim', dest='nextseq_trim',
                        help='nextseq_trim', default=20)
    parser.add_argument(
        '--overlap', help='minimum overlap length, default=5', default=5)
    parser.add_argument('--lowQual', type=int,
                        help='max phred of base as lowQual', default=0)
    parser.add_argument('--lowNum', type=int,
                        help='max number with lowQual allowed', default=2)
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument("--type", help='TCR or BCR', required=True)
    parser.add_argument("--debug", action='store_true')
    parser.add_argument(
        '--iUMI', help='minimum number of UMI of identical receptor type and CDR3')
    parser.add_argument('--rm_files', action='store_true',
                        help='remove redundant fq.gz and bam after running')
    args = vars(parser.parse_args())

    fq_dict, match_dict = parse_map_col4(args['mapfile'], None)

    raw_dir = args['outdir'] + '/data_give/rawdata'
    os.system('mkdir -p %s' % (raw_dir))
    with open(raw_dir + '/ln.sh', 'w') as fh:
        fh.write('cd %s\n' % (raw_dir))
        for s, arr in fq_dict.items():
            fh.write('ln -sf %s %s\n' % (arr[0], s + '_1.fq.gz'))
            fh.write('ln -sf %s %s\n' % (arr[1], s + '_2.fq.gz'))

    logdir = args['outdir'] + '/log'
    os.system('mkdir -p %s' % (logdir))
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    shell = ''
    app = 'celescope'
    mod = args['mod']
    thread = args['thread']
    chemistry = args['chemistry']
    pattern = args['pattern']
    whitelist = args['whitelist']
    linker = args['linker']
    lowQual = args['lowQual']
    lowNum = args['lowNum']
    basedir = args['outdir']
    type = args['type']
    iUMI = args['iUMI']

    assay = __ASSAY__
    steps = __STEPS__
    conda = __CONDA__

    for sample in fq_dict:
        outdir_dic = {}
        index = 0
        for step in steps:
            outdir = f"{basedir}/{sample}/{index:02d}.{step}"
            outdir_dic.update({step: outdir})
            index += 1

        step = "sample"
        cmd = (
            f'{app} {assay} {step} '
            f'--chemistry {chemistry} '
            f'--sample {sample} --outdir {outdir_dic[step]} --assay {assay} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda)
        shell += cmd + '\n'
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
        shell += cmd + '\n'
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
        shell += cmd + '\n'
        last_step = step

        # mapping_vdj
        step = 'mapping_vdj'
        fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--fq {fq} '
            f'--sample {sample} '
            f'--type {type} '
            f'--thread {thread} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--thread {thread} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=15, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell += cmd + '\n'
        last_step = step

        # count_vdj
        step = 'count_vdj'
        UMI_count_filter1_file = f'{outdir_dic["mapping_vdj"]}/{sample}_UMI_count_filtered1.tsv'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} '
            f'--type {type} '
            f'--iUMI {iUMI} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--UMI_count_filter1_file {UMI_count_filter1_file} '
            f'--match_dir {match_dict[sample]} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=8, x=thread)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell += cmd + '\n'
        last_step = step

    # merged report
    step = 'merge_report'
    if mod == 'sjm':
        # add type to steps mapping and count
        for i in range(3, len(steps)):
            steps[i] = f'{type}_{steps[i]}'
        merge_report(
            fq_dict, steps, last_step, sjm_cmd, sjm_order,
            logdir, conda, args['outdir'], args['rm_files']
        )
    if mod == 'shell':
        os.system('mkdir -p ./shell/')
        with open(f'./shell/{sample}.sh', 'w') as f:
            f.write(shell)


if __name__ == '__main__':
    main()
