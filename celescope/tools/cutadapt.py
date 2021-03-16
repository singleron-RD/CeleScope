#!/bin/env python
# coding=utf8

import os
import sys
import re
import subprocess
import logging
from itertools import islice
import pandas as pd
import pysam

from celescope.tools.utils import format_number, log
from celescope.tools.report import reporter
from celescope.tools.Fastq import Fastq


def format_stat(cutadapt_log, samplename):
    fh = open(cutadapt_log, 'r')
    stat_file = os.path.dirname(cutadapt_log) + '/stat.txt'
    # Total reads processed:...Total written (filtered):
    content = islice(fh, 9, 16)
    p_list = []
    for line in content:
        if line.strip() == '':
            continue
        line = re.sub(r'\s{2,}', r'', line)
        line = re.sub(r' bp', r'', line)
        line = re.sub(r'(?<=\d)\s+\(', r'(', line)
        line = line.strip()
        attr = line.split(":")
        p_list.append({"item": attr[0], "value": attr[1]})
    p_df = pd.DataFrame(p_list)
    p_df.iloc[0, 0] = 'Reads with Adapters'
    p_df.iloc[1, 0] = 'Reads too Short'
    p_df.iloc[2, 0] = 'Reads Written'
    p_df.iloc[3, 0] = 'Base Pairs Processed'
    p_df.iloc[4, 0] = 'Base Pairs Quality-Trimmed'
    p_df.iloc[5, 0] = 'Base Pairs Written'
    p_df.to_csv(stat_file, sep=':', index=False, header=None)

    fh.close()


@log
def read_adapter_fasta(adapter_fasta):
    '''
    return ['adapter1=AAA','adapter2=BBB']
    '''
    adapter_args = []
    if adapter_fasta and adapter_fasta!='None':
        with pysam.FastxFile(adapter_fasta) as fh:
            for read in fh:
                adapter_args.append(f'{read.name}={read.sequence}')
    return adapter_args


@log
def cutadapt(args):
    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % (args.outdir))

    adapter_args = read_adapter_fasta(args.adapter_fasta)
    adapter_args += args.adapt

    # run cutadapt
    adapt = []
    for a in adapter_args:
        adapt.append('-a')
        adapt.append(a)

    if not args.not_gzip:
        suffix = ".gz"
    else:
        suffix = ""
    out_fq2 = f'{args.outdir}/{args.sample}_clean_2.fq{suffix}'
    cmd = ['cutadapt'] + adapt + ['-n',
                                  str(len(adapter_args)),
                                  '-j',
                                  str(args.thread),
                                  '-m',
                                  str(args.minimum_length),
                                  '--nextseq-trim=' + str(args.nextseq_trim),
                                  '--overlap',
                                  str(args.overlap),
                                  '-l',
                                  str(args.insert),
                                  '-o',
                                  out_fq2,
                                  args.fq]
    cutadapt.logger.info('%s' % (' '.join(cmd)))
    res = subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    with open(args.outdir + '/cutadapt.log', 'wb') as fh:
        fh.write(res.stdout)

    format_stat(args.outdir + '/cutadapt.log', args.sample)

    t = reporter(
        name='cutadapt',
        assay=args.assay,
        sample=args.sample,
        stat_file=args.outdir +
        '/stat.txt',
        outdir=args.outdir +
        '/..')
    t.get_report()


def get_opts_cutadapt(parser, sub_program):
    if sub_program:
        parser.add_argument('--fq', help='fq file', required=True)
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument('--not_gzip', help="output fastq without gzip", action='store_true')
    parser.add_argument('--adapt',action='append',default=[
            'polyT=A{18}',
            'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',])
    parser.add_argument('--adapter_fasta', help='addtional adapter fasta file')
    parser.add_argument('--minimum_length',dest='minimum_length',help='minimum_length, default=20',default=20)
    parser.add_argument('--nextseq-trim',dest='nextseq_trim',help='nextseq_trim, default=20',default=20)
    parser.add_argument('--overlap',help='minimum overlap length',default=10)
    parser.add_argument('--thread', default=2)
    parser.add_argument('--insert', help="read2 insert length", default=150)


