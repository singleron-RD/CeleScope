#!/bin/env python
#coding=utf8

import argparse
import sys
import os 
import logging
from celescope.tools.__version__ import __VERSION__
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


if __name__ == '__main__':
    main()

def main():
    parser = argparse.ArgumentParser(description='CeleScope')
    parser.add_argument('-v', '--version', action='version', version= __VERSION__)
    subparsers = parser.add_subparsers()

    from celescope.tools.sample_info import sample_info, get_opts_sample
    parser_sample = subparsers.add_parser('sample', description='sample infomation')
    get_opts_sample(parser_sample, True)
    parser_sample.set_defaults(func=sample_info)

    from celescope.tools.barcode import barcode, get_opts_barcode
    parser_barcode = subparsers.add_parser('barcode', description='extract barcode and umi')
    get_opts_barcode(parser_barcode, True)
    parser_barcode.set_defaults(func=barcode)

    from celescope.tools.cutadapt import cutadapt, get_opts_cutadapt
    parser_cutadapt = subparsers.add_parser('cutadapt', description='cutadapt')
    get_opts_cutadapt(parser_cutadapt, True)
    parser_cutadapt.set_defaults(func=cutadapt)

    from celescope.tools.STAR import STAR, get_opts_STAR
    parser_STAR = subparsers.add_parser('STAR')
    get_opts_STAR(parser_STAR, True)
    parser_STAR.set_defaults(func=STAR)

    from celescope.rna_virus.STAR_virus import STAR_virus, get_opts_STAR_virus
    parser_STAR_virus = subparsers.add_parser('STAR_virus')
    get_opts_STAR_virus(parser_STAR_virus, True)
    parser_STAR_virus.set_defaults(func=STAR_virus)

    from celescope.tools.featureCounts import featureCounts, get_opts_featureCounts
    parser_featureCounts = subparsers.add_parser('featureCounts')
    get_opts_featureCounts(parser_featureCounts, True)
    parser_featureCounts.set_defaults(func=featureCounts)

    from celescope.tools.count import count, get_opts_count
    parser_count = subparsers.add_parser('count')
    get_opts_count(parser_count, True)
    parser_count.set_defaults(func=count)

    from celescope.rna_virus.count_virus import count_virus, get_opts_count_virus
    parser_count_virus = subparsers.add_parser('count_virus')
    get_opts_count_virus(parser_count_virus, True)
    parser_count_virus.set_defaults(func=count_virus)   

    from celescope.capture_virus.count_capture_virus import count_capture_virus, get_opts_count_capture_virus
    parser_count_capture_virus = subparsers.add_parser('count_capture_virus')
    get_opts_count_capture_virus(parser_count_capture_virus, True)
    parser_count_capture_virus.set_defaults(func=count_capture_virus)

    from celescope.tools.analysis import analysis, get_opts_analysis
    parser_analysis = subparsers.add_parser('analysis')
    get_opts_analysis(parser_analysis, True)
    parser_analysis.set_defaults(func=analysis) 

    from celescope.rna_virus.analysis_rna_virus import analysis_rna_virus, get_opts_analysis_rna_virus
    parser_analysis_rna_virus = subparsers.add_parser('analysis_rna_virus')
    get_opts_analysis_rna_virus(parser_analysis_rna_virus, True)
    parser_analysis_rna_virus.set_defaults(func=analysis_rna_virus)

    from celescope.capture_virus.analysis_capture_virus import analysis_capture_virus, get_opts_analysis_capture_virus
    parser_analysis_capture_virus = subparsers.add_parser('analysis_capture_virus')
    get_opts_analysis_capture_virus(parser_analysis_capture_virus, True)
    parser_analysis_capture_virus.set_defaults(func=analysis_capture_virus)

    from celescope.rna.rna import rna
    parser_run = subparsers.add_parser('rna', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR(parser_run, False)
    get_opts_featureCounts(parser_run, False)
    get_opts_count(parser_run, False)
    get_opts_analysis(parser_run, False)
    parser_run.set_defaults(func=rna)

    from celescope.rna_virus.rna_virus import rna_virus
    parser_run = subparsers.add_parser('rna_virus', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR(parser_run, False)
    get_opts_STAR_virus(parser_run, False)
    get_opts_featureCounts(parser_run, False)
    get_opts_count_virus(parser_run, False)
    get_opts_analysis_rna_virus(parser_run, False)
    parser_run.set_defaults(func=rna_virus)

    from celescope.capture_virus.capture_virus import capture_virus
    parser_run = subparsers.add_parser('capture_virus', conflict_handler='resolve')
    get_opts_sample(parser_run, False)
    get_opts_barcode(parser_run, False)
    get_opts_cutadapt(parser_run, False)
    get_opts_STAR_virus(parser_run, False)
    get_opts_count_capture_virus(parser_run, False)
    parser_run.set_defaults(func=capture_virus)

    args = parser.parse_args()
    args.func(args)
