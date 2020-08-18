#!/bin/env python
#coding=utf8

import os
import sys
import json
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser('merge report')
parser.add_argument('--samples', help='samples, seperated by comma', required=True)
parser.add_argument('--steps', help='steps', required=True)
args = vars(parser.parse_args())

steps = args['steps'].split(",")
samples = args['samples'].split(',')
result_dict = defaultdict(list)


def pie_label(values, keys):
    total = float(sum(values))
    return ['%s:%.2f%%'%(k.replace(' Regions',''), v/total*100) for k, v in zip(keys, values)]


summarys = [step + '_summary' for step in steps]
for sample in samples:
    data_json = f'./{sample}/.data.json'
    data_dic = json.load(open(data_json))
    for summary in summarys:
        if summary not in data_dic.keys():
            continue
        # add title
        if sample == samples[0]:
            result_dict[summary].append('\t'.join([x[0].replace(' ', '_') for x in data_dic[summary]]))
        result_dict[summary].append('\t'.join([x[1].replace(' ', '') for x in data_dic[summary]]))
    
with open('./merge.xls', 'w') as fh:
    for summary in summarys:
        if summary not in data_dic.keys():
            continue
        fh.write('##' + summary + '\n')
        for k in result_dict[summary]:
            fh.write(k + '\n')
        fh.write('\n')

