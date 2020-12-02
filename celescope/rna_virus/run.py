#!/bin/env python
# coding=utf8
import os
import re
import io
import logging
import gzip
import json
import subprocess
from collections import defaultdict, namedtuple
import sys


def run(args):
    steps = ['sample', 'barcode', 'cutadapt', 'STAR', "STAR_virus", "featureCounts", "count", "count_virus",
             'analysis']
    sample = args.sample
    args.assay = "rna_virus"

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

    step = "STAR"
    args.fq = f'{outdir_dic["cutadapt"]}/{sample}_clean_2.fq.gz'
    args.outdir = f'{outdir_dic[step]}/'
    args.out_unmapped = True
    from celescope.tools.STAR import STAR
    STAR(args)

    step = "STAR_virus"
    args.input_read = f'{outdir_dic["STAR"]}/{sample}_Unmapped.out.mate1'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.rna_virus.STAR_virus import STAR_virus
    STAR_virus(args)

    step = 'featureCounts'
    args.input = f'{outdir_dic["STAR"]}/{sample}_Aligned.sortedByCoord.out.bam'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.tools.featureCounts import featureCounts
    featureCounts(args)

    step = 'count'
    args.bam = f'{outdir_dic["featureCounts"]}/{sample}_name_sorted.bam'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.tools.count import count
    count(args)

    step = 'count_virus'
    args.virus_bam = f'{outdir_dic["STAR_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
    args.barcode_file = f'{outdir_dic["count"]}/matrix_10X/{sample}_cellbarcode.tsv'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.virus.count_virus import count_virus
    count_virus(args)

    step = 'analysis'
    args.matrix_file = f'{outdir_dic["count"]}/{sample}_matrix.xls'
    args.virus_file = f'{outdir_dic["virus_count"]}/{sample}_virus_UMI_count.tsv'
    args.outdir = f'{outdir_dic[step]}/'
    from celescope.tools.analysis import analysis
    analysis(args)
