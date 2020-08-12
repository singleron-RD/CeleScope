#!/bin/env python
#coding=utf8
import os, re, io, logging, gzip, json
import subprocess
from collections import defaultdict, namedtuple
import sys


def virus(args):
    sample = args.sample
    baseDir = sample

    if not os.path.exists(baseDir):
        os.system('mkdir -p %s' % baseDir)   

    args.outdir = baseDir + '/00.sample'    
    from tools.sample_info import sample_info
    sample_info(args)    

    args.outdir = baseDir + '/01.barcode'    
    from tools.barcode import barcode
    args.fq = baseDir + '/01.barcode/' + sample + '_2.fq.gz'
    barcode(args)

    args.outdir = baseDir + '/02.cutadapt'
    from tools.cutadapt import cutadapt
    cutadapt(args)

    args.thread = 4
    args.clean_read = baseDir + '/02.cutadapt/' + sample + '_clean_2.fq.gz'
    args.outdir = baseDir + '/03.STAR'
    from tools.STAR import STAR
    STAR(args)

    args.clean_read = baseDir + '/03.STAR/' + sample + '_Unmapped.out.mate1'
    args.outdir = baseDir + '/03.STAR'
    from virus.STAR_virus import STAR_virus
    STAR_virus(args) 

    args.input = baseDir + '/03.STAR/' + sample + '_Aligned.sortedByCoord.out.bam'
    args.outdir = baseDir + '/04.featureCounts'
    from tools.featureCounts import featureCounts
    featureCounts(args)

    args.bam = baseDir + '/04.featureCounts/' + sample + '_name_sorted.bam'
    args.virus_bam = baseDir + '/03.STAR/' + sample + '_virus_Aligned.sortedByCoord.out.bam'
    args.outdir = baseDir + '/05.count'
    from virus.count_virus import count_virus
    count_virus(args)

    args.matrix_file = baseDir + '/05.count/' + sample + '_matrix.xls'
    args.virus_file = baseDir + '/05.count/' + sample + '_virus_UMI_count.tsv'
    args.outdir = baseDir + '/06.analysis'
    from tools.analysis import analysis
    analysis(args)   


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Single cell RNA-Seq virus detection')
    args = parser.parse_args()
    virus(args)


if __name__ == '__main__':
    main()
