# !/bin/env python
# coding=utf8
from celescope.smk.__init__ import __STEPS__, __ASSAY__


def run(args):
    steps = __STEPS__
    args.assay = __ASSAY__
    sample = args.sample

    outdir_dic = {}
    index = 0
    for step in steps:
        outdir = f"{sample}/{index:02d}.{step}"
        outdir_dic.update({step: outdir})
        index += 1

    step = "sample"
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.tools.sample_info import sample_info
    sample_info(args)

    step = "barcode"
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.tools.barcode import barcode
    barcode(args)

    step = "cutadapt"
    args.outdir = f'{outdir_dic[step]}/'
    args.fq = f'{outdir_dic["barcode"]}/{sample}_2.fq.gz'
    from celescope.tools.cutadapt import cutadapt
    cutadapt(args)

    step = 'demultiplex'
    args.outdir = f'{outdir_dic[step]}/'
    args.SMK_read2 = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
    from celescope.smk.demultiplex import demultiplex
    demultiplex(args)

