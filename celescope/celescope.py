#!/bin/env python
#coding=utf8

import argparse
import sys
import os 
"""
toolsdir = os.path.realpath(sys.path[0] + '/tools')
rnadir = os.path.realpath(sys.path[0] + '/rna')
virusdir = os.path.realpath(sys.path[0] + '/virus')
sys.path.append(toolsdir)
sys.path.append(rnadir)
#sys.path.append(virusdir)
"""
from tools.__version__ import __VERSION__


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CeleScope')
    parser.add_argument('-v', '--version', action='version', version= __VERSION__)
    subparsers = parser.add_subparsers()

    from tools.sample_info import sample_info, get_opts_sample
    parser_sample = subparsers.add_parser('sample', description='sample infomation')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    from tools.barcode import barcode, get_opts_barcode
    parser_barcode = subparsers.add_parser('barcode', description='extract barcode and umi')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    from tools.cutadapt import cutadapt, get_opts_cutadapt
    parser_cutadapt = subparsers.add_parser('cutadapt', description='cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from tools.STAR import STAR, get_opts_STAR
    parser_STAR = subparsers.add_parser('STAR')
    get_opts_STAR(parser_STAR, True)
    parser_STAR.set_defaults(func=STAR)

    from virus.STAR_virus import STAR_virus, get_opts_STAR_virus
    parser_STAR_virus = subparsers.add_parser('STAR_virus')
    get_opts_STAR_virus(parser_STAR_virus, True)
    parser_STAR_virus.set_defaults(func=STAR_virus)

    from tools.featureCounts import featureCounts, get_opts_featureCounts
    parser_featureCounts = subparsers.add_parser('featureCounts')
    get_opts_featureCounts(parser_featureCounts, True)
    parser_featureCounts.set_defaults(func=featureCounts)

    from tools.count import count, get_opts_count
    parser_count = subparsers.add_parser('count')
    get_opts_count(parser_count, True)
    parser_count.set_defaults(func=count)

    from virus.count_virus import count_virus, get_opts_count_virus
    parser_count_virus = subparsers.add_parser('count_virus')
    get_opts_count_virus(parser_count_virus, True)
    parser_count_virus.set_defaults(func=count_virus)    

    from tools.analysis import analysis, get_opts_analysis
    parser_analysis = subparsers.add_parser('analysis')
    get_opts_analysis(parser_analysis, True)
    parser_analysis.set_defaults(func=analysis)

    from rna.rna import rna
    parser_run = subparsers.add_parser('rna', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR(parser_run, False)
    get_opts_featureCounts(parser_run, False)
    get_opts_count(parser_run, False)
    get_opts_analysis(parser_run, False)
    parser_run.set_defaults(func=rna)

    from virus.virus import virus
    parser_run = subparsers.add_parser('virus', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR(parser_run, False)
    get_opts_STAR_virus(parser_run, False)
    get_opts_featureCounts(parser_run, False)
    get_opts_count_virus(parser_run, False)
    get_opts_analysis(parser_run, False)
    parser_run.set_defaults(func=virus)

    args = parser.parse_args()
    args.func(args)
