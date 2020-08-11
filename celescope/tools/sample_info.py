#!/bin/env python
# coding=utf8

import os
import pandas as pd
import sys
from tools.utils import getlogger
from tools.__version__ import __VERSION__
from tools.report import reporter

logger1 = getlogger()

def sample_info(args):
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % args.outdir)
    sample = args.sample
    description = args.description
    version = args.version
    outdir = args.outdir
    chemistry = args.chemistry
    if not chemistry:
        chemistry = "Customized"
    #transcriptome = args.genomeDir.split("/")[-1]

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    stat = pd.DataFrame({"item": ["Sample ID", "Description", "Chemistry", "Software Version"],
        "count": [sample, description, chemistry, version]}, columns=["item", "count"])
    stat_file = outdir + "/stat.txt"
    stat.to_csv(stat_file, sep=":", header=None, index=False)

    t = reporter(name='sample', sample=args.sample, stat_file=stat_file, outdir=outdir + '/..')
    t.get_report()
    logger1.info("Generating sample info done.")


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--genomeDir', help='genomeDir', required=True)
    parser.add_argument('--description', help='sample description', default="scRNA-Seq")
    parser.add_argument('--chemistry', help='chemistry')
    parser.add_argument('--version', help='software version', default=__VERSION__)