#!/bin/env python
#coding=utf8


def run(args):
    steps = ['sample', 'barcode', 'cutadapt',  "STAR_fusion", "count_fusion"]
    sample = args.sample
    args.assay = "fusion"

    outdir_dic = {}
    index = 0
    for step in steps:
        outdir = f"{sample}/{index:02d}.{step}"
        outdir_dic.update({step: outdir})
        index += 1

    step = "sample"
    args.outdir = f'{sample}/{outdir_dic["step"]}/'
    from celescope.tools.sample_info import sample_info
    sample_info(args)   

    step = "barcode"
    args.outdir = f'{sample}/{outdir_dic["step"]}/'   
    from celescope.tools.barcode import barcode
    barcode(args)

    step = "cutadapt"
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    args.fq = f'{outdir_dic["barcode"]}/{sample}_2.fq.gz'
    from celescope.tools.cutadapt import cutadapt
    cutadapt(args)

    step = "STAR_fusion"
    args.input_read = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from celescope.fusion.STAR_fusion import STAR_fusion
    STAR_fusion(args)  

    step = 'count_fusion'
    args.virus_bam = f'{outdir_dic["STAR_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from celescope.fusion.count_fusion import count_fusion
    count_fusion(args)