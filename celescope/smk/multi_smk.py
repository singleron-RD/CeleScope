import os
import glob
import sys
import argparse
import re
from collections import defaultdict
from celescope.__init__ import __CONDA__
from celescope.smk.__init__ import __STEPS__, __ASSAY__
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
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument(
        "--UMI_min",
        help="cells have SMK_UMI>=UMI_min are considered as valid cell",
        default="auto")
    parser.add_argument("--dim", help="SMK tag dimension", default=1)
    parser.add_argument(
        "--SNR_min",
        help="minimum signal to noise ratio",
        default="auto")
    parser.add_argument("--SMK_pattern", help="SMK read2 pattern")
    parser.add_argument("--SMK_linker", help="SMK read2 linker fasta path")
    parser.add_argument("--SMK_barcode", help="SMK read2 barcode fasta path ")
    args = parser.parse_args()

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
    fq_dict, match_dict = parse_map_col4(args.mapfile, "auto")

    # link
    link_data(outdir, fq_dict)

    # custom args
    thread = args.thread
    UMI_min = args.UMI_min
    dim = args.dim
    SNR_min = args.SNR_min
    SMK_pattern = args.SMK_pattern
    SMK_linker = args.SMK_linker
    SMK_barcode = args.SMK_barcode

    # mk log dir
    logdir = outdir + '/log'
    os.system('mkdir -p %s' % (logdir))

    # script init
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    shell = ''

    # run
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

        # mapping_smk
        step = 'mapping_smk'
        SMK_read2 = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--SMK_read2 {SMK_read2} '
            f'--SMK_pattern {SMK_pattern} '
            f'--SMK_barcode {SMK_barcode} '
            f'--SMK_linker {SMK_linker} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell += cmd + '\n'
        last_step = step

        # count_smk
        step = 'count_smk'
        read_file = f'{outdir_dic["mapping_smk"]}/{sample}_read_count.tsv'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--match_dir {match_dict[sample]} '
            f'--read_file {read_file} '
            f'--dim {dim} '
            f'--UMI_min {UMI_min} '
            f'--SNR_min {SNR_min} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell += cmd + '\n'
        last_step = step

        # analysis_smk
        step = 'analysis_smk'
        tsne_tag_file = f'{outdir_dic["count_smk"]}/{sample}_tsne_tag.tsv'
        cmd = (
            f'{app} {assay} {step} '
            f'--sample {sample} '
            f'--outdir {outdir_dic[step]} '
            f'--assay {assay} '
            f'--match_dir {match_dict[sample]} '
            f'--tsne_tag_file {tsne_tag_file} '
        )
        sjm_cmd += generate_sjm(cmd, f'{step}_{sample}', conda, m=5, x=1)
        sjm_order += f'order {step}_{sample} after {last_step}_{sample}\n'
        shell += cmd + '\n'
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
        with open(f'./shell/{sample}.sh', 'w') as f:
            f.write(shell)


if __name__ == '__main__':
    main()
