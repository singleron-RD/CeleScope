#!/bin/env python
#coding=utf8
import os, re, io, logging, gzip, json
import subprocess
from collections import defaultdict, namedtuple
import sys


def capture_virus(args):
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

    args.input_read = baseDir + '/02.cutadapt/' + sample + '_clean_2.fq.gz'
    args.outdir = baseDir + '/03.STAR_virus'
    from virus.STAR_virus import STAR_virus
    STAR_virus(args) 

    args.virus_bam = baseDir + '/03.STAR_virus/' + sample + '_virus_Aligned.sortedByCoord.out.bam'
    args.outdir = baseDir + '/04.count'
    from virus.capture_count_virus import capture_count_virus
    capture_count_virus(args)
 