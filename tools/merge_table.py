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
parser.add_argument('--workdir', help='working dir', required=True)
args = vars(parser.parse_args())

samples = args['samples'].split(',')
data_json = [args['workdir'] + '/' + s + '/' + '.data.json' for s in samples]
result_dict = defaultdict(list)

def pie_label(values, keys):
    total = float(sum(values))
    return ['%s:%.2f%%'%(k.replace(' Regions',''), v/total*100) for k, v in zip(keys, values)]

for n, i in zip(samples, data_json):
    tmp = json.load(open(i))
    for j in ['sample_summary','barcode_summary', 'cutadapt_summary', 'STAR_summary', 'featureCounts_summary', 'count_summary']:
        if n==samples[0]: 
            result_dict[j].append('\t'.join([x[0].replace(' ','') for x in tmp[j]]))
        result_dict[j].append('\t'.join([x[1].replace(' ','') for x in tmp[j]]))
    
with open(args['workdir']+'/merge.xls', 'w') as fh:
    for j in ['sample_summary','barcode_summary', 'cutadapt_summary', 'STAR_summary', 'featureCounts_summary', 'count_summary']:
        fh.write('##' + j + '\n')
        for k in result_dict[j]:
            fh.write(k + '\n')
        fh.write('\n')

