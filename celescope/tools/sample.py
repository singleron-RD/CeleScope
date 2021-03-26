#!/bin/env python
# coding=utf8

import os
import pandas as pd
import sys
from celescope.__init__ import __VERSION__, ASSAY_DICT
from celescope.tools.utils import *
from celescope.tools.report import reporter
from celescope.tools.__init__ import __PATTERN_DICT__
from .Chemistry import Chemistry


@add_log
def sample(args):

    sample_name = args.sample
    assay = args.assay
    assay_description = ASSAY_DICT[assay]
    version = __VERSION__
    outdir = args.outdir
    chemistry = args.chemistry

    # get chemistry
    if chemistry == 'auto':
        fq1 = args.fq1
        ch = Chemistry(fq1)
        chemistry = ch.check_chemistry()
        chemistry = ",".join(set(chemistry))
    else:
        chemistry = args.chemistry
    

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    stat = pd.DataFrame({
        "item": ["Sample ID", "Assay", "Chemistry", "Software Version"],
        "count": [sample_name, assay_description, chemistry, version],
        },
        columns=["item", "count"]
    )
    stat_file = outdir + "/stat.txt"
    stat.to_csv(stat_file, sep=":", header=None, index=False)

    t = reporter(
        name='sample',
        assay=assay,
        sample=sample_name,
        stat_file=stat_file,
        outdir=outdir + '/..')
    t.get_report()
    return chemistry


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='read1 fq file')
    parser.add_argument('--chemistry', choices=__PATTERN_DICT__.keys(), help='chemistry version', default='auto')
    return parser
    
