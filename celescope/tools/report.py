#!/bin/env python
#coding=utf8

import os, sys, json
import argparse
import pandas as pd
import io
from jinja2 import Environment, PackageLoader, select_autoescape, FileSystemLoader

env = Environment(
    loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
    autoescape=select_autoescape(['html', 'xml'])
)


class reporter:
    def __init__(self, assay, name, outdir, sample, stat_file=None, plot=None):
        self.assay = assay
        self.name = name
        self.stat_file = stat_file
        self.outdir = outdir
        self.sample = sample
        self.plot = plot
  
    def get_report(self):
        template = env.get_template(f'html/{self.assay}/base.html')
        json_file = self.outdir + '/.data.json'
        if not os.path.exists(json_file):
            data = {}
        else:
            fh = open(json_file)
            data = json.load(fh)
            #if 'STAR_plot' in data:
                #data['STAR_plot']['region_labels'] = json.dumps(data['STAR_plot']['region_labels'])
            fh.close()

        if self.stat_file:
            df = pd.read_table(self.stat_file, header=None, sep=':', dtype=str)
            #data[self.name + '_summary'] = df.T.values.tolist()
            data[self.name + '_summary'] = df.values.tolist()
  
        if self.plot:
            data[self.name + '_plot'] = self.plot

        report_html = "{outdir}/{sample}_report.html".format(outdir=self.outdir,sample=self.sample)
        with io.open(report_html, 'w',encoding='utf8') as fh:
            html = template.render(data)
            #fh.write(html.encode('utf-8'))
            fh.write(html)

        with open(json_file, 'w') as fh:
            json.dump(data, fh)






