#!/bin/env python
#coding=utf8
import os, re, io, logging, gzip, json
import subprocess
from collections import defaultdict, namedtuple
import sys
sys.path.append('../tools')


def rna(args):
    #tmp = vars(args)
    sample = args.sample
    baseDir = args.sample  
    args.assay = "rna"  

    args.outdir = baseDir + '/00.sample'
    args.description = "Single Cell RNA-Seq"
    from celescope.tools.sample_info import sample_info
    sample_info(args)

    args.outdir = baseDir + '/01.barcode'    
    from celescope.tools.barcode import barcode
    barcode(args)

    args.fq = baseDir + '/01.barcode/' + sample + '_2.fq.gz'
    args.outdir = baseDir + '/02.cutadapt'
    from celescope.tools.cutadapt import cutadapt
    cutadapt(args)

    args.fq = baseDir + '/02.cutadapt/' + sample + '_clean_2.fq.gz'
    args.outdir = baseDir + '/03.STAR'
    from celescope.tools.STAR import STAR
    STAR(args)

    args.input = baseDir + '/03.STAR/' + sample + '_Aligned.sortedByCoord.out.bam'
    args.outdir = baseDir + '/04.featureCounts'
    from celescope.tools.featureCounts import featureCounts
    featureCounts(args)

    args.bam = baseDir + '/04.featureCounts/' + sample + '_name_sorted.bam'
    args.outdir = baseDir + '/05.count'
    from celescope.tools.count import count
    count(args)

    args.matrix_file = baseDir + '/05.count/' + sample + '_matrix.xls'
    args.outdir = baseDir + '/06.analysis'
    from celescope.tools.analysis import analysis
    analysis(args)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Single cell RNA-Seq')
    args = parser.parse_args()
    rna(args)


if __name__ == '__main__':
    main()
