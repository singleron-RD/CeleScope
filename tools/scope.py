#!/bin/env python
#coding=utf8

import argparse
import sys


__VERSION__ = "CeleScope v1.1.0"


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CeleScope')
    parser.add_argument('-v', '--version', action='version', version= __VERSION__)
    subparsers = parser.add_subparsers()

    from sample_info import sample_info, get_opts_sample
    parser_sample = subparsers.add_parser('sample', description='sample infomation')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    from barcode import barcode, get_opts_barcode
    parser_barcode = subparsers.add_parser('barcode', description='extract barcode and umi')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    from cutadapt import cutadapt, get_opts_cutadapt
    parser_cutadapt = subparsers.add_parser('cutadapt', description='cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from STAR import STAR, get_opts_STAR
    parser_STAR = subparsers.add_parser('STAR')
    get_opts_STAR(parser_STAR, True)
    parser_STAR.set_defaults(func=STAR)

    from featureCounts import featureCounts, get_opts_featureCounts
    parser_featureCounts = subparsers.add_parser('featureCounts')
    get_opts_featureCounts(parser_featureCounts, True)
    parser_featureCounts.set_defaults(func=featureCounts)

    from count import count, get_opts_count
    parser_count = subparsers.add_parser('count')
    get_opts_count(parser_count, True)
    parser_count.set_defaults(func=count)

    from analysis import analysis, get_opts_analysis
    parser_analysis = subparsers.add_parser('analysis')
    get_opts_analysis(parser_analysis, True)
    parser_analysis.set_defaults(func=analysis)

    from rna import rna
    parser_run = subparsers.add_parser('RNA', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR(parser_run, False)
    get_opts_featureCounts(parser_run, False)
    get_opts_count(parser_run, False)
    get_opts_analysis(parser_run, False)
    parser_run.set_defaults(func=rna)

    args = parser.parse_args()
    args.func(args)
