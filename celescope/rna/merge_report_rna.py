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
    
    # pie
    patches, texts = plt.pie(tmp['STAR_plot']['region_values'])
    plt.legend(patches, pie_label(tmp['STAR_plot']['region_values'], tmp['STAR_plot']['region_labels']),bbox_to_anchor=(1,0.5), loc="center right")
    plt.savefig(n+'_pie.png', bbox_inches="tight", dpi=300)
    plt.clf()

    #ax.pie(tmp['STAR_plot']['region_values'], labels=tmp['STAR_plot']['region_labels'])

with open(args['workdir']+'/merge.xls', 'w') as fh:
    for j in ['barcode_summary', 'cutadapt_summary', 'STAR_summary', 'featureCounts_summary', 'count_summary']:
        fh.write('##' + j + '\n')
        for k in result_dict[j]:
            fh.write(k + '\n')
        fh.write('\n')

'''
plt.pie(x['STAR_plot']['region_values'], labels=x['STAR_plot']['region_labels'])

def str2list(s): 
    return [int(i) for i in s.replace('[','').replace(']','').split(',')]

cells = str2list(x['Cells'])
Background = str2list(x['Background'])
plt.loglog(list(range(len(cells)+1,len(cells)+1+len(Background))),Background)
plt.savefig('xxx')

plt.plot(x['percentile'], x['Saturation'])
plt.plot(x['percentile'], x['MedianGeneNum'])
'''
