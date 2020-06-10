#!/bin/env python
#coding=utf8

import os, sys, json
import argparse
import pandas as pd
import io
from jinja2 import Environment, PackageLoader, select_autoescape, FileSystemLoader
env = Environment(
    loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates'),
    autoescape=select_autoescape(['html', 'xml'])
)

class reporter:
    def __init__(self, name, outdir, stat_file=None, plot=None):
        self.name = name
        self.stat_file = stat_file
        self.outdir = outdir
        self.plot = plot
  
    def get_report(self):
        template = env.get_template('base.html')
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

        with io.open(self.outdir + '/report.html', 'w',encoding='utf8') as fh:
            html = template.render(data)
            #fh.write(html.encode('utf-8'))
            fh.write(html)

        with open(json_file, 'w') as fh:
            json.dump(data, fh)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='report')
    parser.add_argument('--basedir', help='output base dir',required=True)
    parser.add_argument('--sample', help='sample name', required=True)
    args = parser.parse_args()

    basedir = args.basedir
    sample = args.sample

    indexs = ['00','01','02','03','04','05','06']
    steps = ['sample','barcode','cutadapt','STAR','featureCounts','count','analysis']
    for i in range(len(steps)):

        index = indexs[i]
        step = steps[i]
        outdir = basedir + '/{index}.{step}'.format(index=index,step=step)
        stat_file = outdir + "/stat.txt"
        
        if step=="STAR":
            from STAR import format_stat
            region_txt = outdir + '/' + sample + '_region.log'
            plot = format_stat(outdir+'/'+sample+'_Log.final.out', region_txt, sample)
            t = reporter(name=step, stat_file=stat_file, outdir=outdir + '/..',plot=plot)
            t.get_report()
            continue

        if step == "featureCounts":
            from featureCounts import format_stat
            format_stat(outdir+'/'+sample+'.summary', sample)

        if step == "count":
            from count import  report_prepare
            marked_counts_file = outdir + '/' + sample + '_counts.txt'
            downsample_file = outdir + '/' + sample + '_downsample.txt'
            report_prepare(marked_counts_file, downsample_file, outdir + '/..')

        if step == "analysis":
            from analysis import  report_prepare
            tsne_df_file = "{outdir}/tsne_coord.tsv".format(outdir=outdir)
            marker_df_file = "{outdir}/markers.tsv".format(outdir=outdir)
            tsne_df = pd.read_csv(tsne_df_file,sep="\t")
            marker_df = pd.read_csv(marker_df_file,sep="\t")
            report_prepare(outdir,tsne_df,marker_df)
            t = reporter(name='analysis', outdir=outdir + '/..')
            t.get_report()
            continue

        t = reporter(name=step, stat_file=stat_file, outdir=outdir + '/..')
        t.get_report()




