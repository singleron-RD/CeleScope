#!/bin/env python
#coding=utf8
import os, re, io, logging, gzip, json
import subprocess
from collections import defaultdict, namedtuple
logging.basicConfig(format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

"""
def get_opts6(parser,sub_program): 
    parser.add_argument('--skip', help='step and steps after will not run, eg. STAR,featureCounts,count', default='')
"""

def run(args):
    #tmp = vars(args)
    sample = args.sample
    baseDir = args.sample    

    args.outdir = baseDir + '/00.sample'
    from sample_info import sample_info
    sample_info(args)
 
    args.outdir = baseDir + '/01.barcode'    
    from barcode import barcode
    barcode(args)

    args.fq = baseDir + '/01.barcode/' + sample + '_2.fq.gz'
    args.outdir = baseDir + '/02.cutadapt'
    from cutadapt import cutadapt
    cutadapt(args)

    args.fq = baseDir + '/02.cutadapt/' + sample + '_clean_2.fq.gz'
    args.outdir = baseDir + '/03.STAR'
    args.runThreadN = 6
    from STAR import STAR
    STAR(args)

    args.input = baseDir + '/03.STAR/' + sample + '_Aligned.sortedByCoord.out.bam'
    args.outdir = baseDir + '/04.featureCounts'
    args.runThreadN = 6
    from featureCounts import featureCounts
    featureCounts(args)

    args.bam = baseDir + '/04.featureCounts/' + sample + '_name_sorted.bam'
    args.outdir = baseDir + '/05.count'
    from count import count
    count(args)

    args.matrix_file = baseDir + '/05.count/' + sample + '_matrix.xls'
    args.outdir = baseDir + '/06.analysis'
    from analysis import analysis
    analysis(args)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='run all steps')
    args = parser.parse_args()
    run(args)

if __name__ == '__main__':
    main()
