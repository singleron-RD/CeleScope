#!/bin/env python
# coding=utf8

from collections import defaultdict
import os
import json
import argparse
from celescope.tools.utils import log


@log
def merge_report():
    parser = argparse.ArgumentParser('merge report')
    parser.add_argument('--outdir', help='outdir', required=True)
    parser.add_argument(
        '--samples', help='samples, seperated by comma', required=True)
    parser.add_argument('--steps', help='steps', required=True)
    parser.add_argument('--rm_files', action='store_true',
                        help='remove all fq.gz and bam after running')
    args = vars(parser.parse_args())

    outdir = args['outdir']
    os.chdir(outdir)
    steps = args['steps'].split(",")
    samples = args['samples'].split(',')
    result_dict = defaultdict(list)

    summarys = [step + '_summary' for step in steps]
    for sample in samples:
        data_json = f'./{sample}/.data.json'
        data_dic = json.load(open(data_json))
        for summary in summarys:
            if summary not in data_dic.keys():
                continue
            # add title
            if sample == samples[0]:
                result_dict[summary].append(
                    '\t'.join([str(x[0]).replace(' ', '_') for x in data_dic[summary]])
                )
            result_dict[summary].append(
                '\t'.join([str(x[1]).replace(' ', '') for x in data_dic[summary]])
            )

    with open('./merge.xls', 'w') as fh:
        for summary in summarys:
            if summary not in data_dic.keys():
                continue
            fh.write('##' + summary + '\n')
            for k in result_dict[summary]:
                fh.write(k + '\n')
            fh.write('\n')

    if args['rm_files']:
        rm_files()


@log
def rm_files():
    cmd = '''
        find . -iname '*.fq*' -delete;
        find . -iname '*.bam' -not -path './*/*.featureCounts/*name_sorted.bam' -delete;
    '''
    os.system(cmd)


if __name__ == '__main__':
    merge_report()
