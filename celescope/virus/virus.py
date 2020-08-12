#!/bin/env python
#coding=utf8
import os, re, io, logging, gzip, json
import subprocess
from collections import defaultdict, namedtuple
import sys


def virus(args):
    steps = ['sample', 'barcode', 'cutadapt', 'STAR', "STAR_virus", "featureCounts", "count", "count_virus",
        'analysis']
    sample = args.sample

    if not os.path.exists(baseDir):
        os.system('mkdir -p %s' % baseDir)   

    outdir_dic = {}
    index = 0
    for step in steps:
        outdir = f"{baseDir}/{sample}/{index:02d}.{step}"
        outdir_dic.update({step: outdir})
        index += 1

    step = "sample"
    args.outdir = f'{sample}/{outdir_dic["step"]}/'
    from tools.sample_info import sample_info
    sample_info(args)   

    step = "barcode"
    args.outdir = f'{sample}/{outdir_dic["step"]}/'   
    from tools.barcode import barcode
    barcode(args)

    step = "cutadapt"
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    args.fq = f'{outdir_dic["barcode"]}/{sample}_2.fq.gz'
    from tools.cutadapt import cutadapt
    cutadapt(args)

    step = "STAR"
    args.fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from tools.STAR import STAR
    STAR(args)

    step = "STAR_virus"
    args.input_read = f'{outdir_dic["STAR"]}/{sample}_Unmapped.out.mate1'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from virus.STAR_virus import STAR_virus
    STAR_virus(args) 

    step = 'featureCounts'
    args.input = f'{outdir_dic["STAR"]}/{sample}_Aligned.sortedByCoord.out.bam'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from tools.featureCounts import featureCounts
    featureCounts(args)

    step = 'count'
    args.bam = f'{outdir_dic["featureCounts"]}/{sample}_name_sorted.bam'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from tools.count import count
    count(args)

    step = 'count_virus'
    args.virus_bam = f'{outdir_dic["STAR"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
    args.barcode_file = f'{outdir_dic["STAR"]}/matrix_10X/{sample}_cellbarcode.tsv'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from virus.count_virus import count_virus
    count_virus(args)

    step = 'analysis'
    args.matrix_file = f'{outdir_dic["count"]}/{sample}_matrix.xls'
    args.outdir = f'{sample}/{outdir_dic["step"]}/' 
    from tools.analysis import analysis
    analysis(args)   

