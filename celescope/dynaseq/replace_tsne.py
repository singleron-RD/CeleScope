#!/bin/env python
# coding=utf8

import os
import sys
import json
import logging
import argparse
import re
import numpy as np
import pandas as pd
import glob
import plotly
import plotly.graph_objects as go
from celescope.tools.report import reporter
from celescope.tools.utils import *
from celescope.tools.Step import Step, s_common

toolsdir = os.path.dirname(__file__)


@add_log
def dot_tsne(repfile,tsnefile,outfile):
    cells = {}
    with open(repfile, 'r') as f:
        for i in f:
            ii = i.strip().split()
            cells[ii[0]] = ii[1]

    outf = open(outfile, 'w')
    outf.write("Cell\ttSNE_1\ttSNE_2\tCluster\tratio\n")
    with open(tsnefile, 'r') as f:
        f.readline()
        for i in f:
            ii = i.strip().split()
            if ii[0] in cells:
                outl = '\t'.join(ii[0:4])+'\t'+cells[ii[0]]+'\n'
            else:
                outl = '\t'.join(ii[0:4])+'\t0'+'\n'
            outf.write(outl)
    outf.close()

@add_log
def tsne_plot(txt,outdiv,outhtml):
    import plotly.express as px
    import plotly.graph_objects as go
    import plotly

    df = pd.read_table(txt)
    df.sort_values(by="ratio")
    newtitle="t-SNE plot Colored by RNA Turn-over rate"

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df['tSNE_1'], y=df['tSNE_2'], mode='markers',
                         marker_opacity=0.9,marker_size=4,marker_color=df['ratio'],
                         marker_colorscale="PuBu", marker_showscale=True,
    ))
    fig.update_layout(height=600, width=600,title_text=newtitle)
    fig.update_layout(plot_bgcolor = '#FFFFFF')
    fig.update_xaxes(showgrid=False,linecolor='black', showline=True, ticks='outside',title_text='t-SNE1')
    fig.update_yaxes(showgrid=False,linecolor='black', showline=True, ticks='outside',title_text='t-SNE2')

    div = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
    '''
    with open(outdiv, 'w') as fno:
        fno.write("<h3>RNA Turn-over rate in clusters</h3>\n")
        fno.write(div + '\n')
    fig.write_html(outhtml)
    '''
    return div


def tsne_table(txt):
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix

    marker_gene_table = txt.to_html(
            escape=False,
            index=False,
            table_id='replacement_table_cluster',
            justify="center")

    return marker_gene_table


def file_stat(infile,clu):
    clus = list(set(clu.values()))
    cluster = {}
    for c in clus:
        cluster[c] = {}
    fn = open(infile,"r")
    fnh = fn.readline().strip().split()
    for i in fn:
        ii = i.strip().split()
        for j in range(1,len(ii)):
            if ii[j] == 'NA': continue
            if fnh[j] not in clu: continue
            if ii[0] not in cluster[clu[fnh[j]]]:
                cluster[clu[fnh[j]]][ii[0]] = []
            cluster[clu[fnh[j]]][ii[0]].append(float(ii[j]))
    fn.close()
    return cluster

def tsne_file(infile):
    clu = {}
    with open(infile) as f:
        f.readline()
        for i in f:
            ii = i.strip().split()
            clu[ii[0]] = ii[3]
    return clu


@add_log
def top_gene_cluster(matrix,tsnefile,outfile,mincell=5,topgene=10):
    tsne = tsne_file(tsnefile)
    cluster = file_stat(matrix,tsne)

    w = open(outfile,'w')
    w.write("cluster\tgene\treplacement_rate\tcells\n")
    for c in cluster:
        tmp = {}
        for g in cluster[c]:
            gt = sum(cluster[c][g]) / len(cluster[c][g])
            tmp[g] = gt
        sorttmp = sorted(tmp.items(), key=lambda item:item[1], reverse=True)
        tmpn = 0
        for x in sorttmp:
            if len(cluster[c][x[0]]) < mincell: continue
            tmpn += 1
            if tmpn > topgene: break
            w.write('cluster'+c+'\t'+x[0]+'\t'+str(x[1])+'\t'+str(len(cluster[c][x[0]]))+'\n')
    w.close()


def report_prepare(outdiv, outable, step_obj):
    step_obj.add_data_item(replace_tsne=outdiv)
    step_obj.add_data_item(replace_tsne_table=outable)



@add_log
def replace_tsne(args):

    step_name = "replace_tsne"
    step = Step(args, step_name)

    # check dir
    outdir = args.outdir
    sample = args.sample
    assay = args.assay
    tsnefile = args.tsne 
    matfile = args.mat
    repfile = args.rep

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))
    outdot = os.path.join(outdir, sample+'.rep_in_tsne.txt')
    outhtml = os.path.join(outdir, sample+'.tsne.html')
    outdiv = os.path.join(outdir, sample+'.tsne.div')
    outtbl = os.path.join(outdir, sample+'.rep_in_tsne_top10.txt')


    # rep in cells in cluster
    dot_tsne(repfile,tsnefile,outdot)
    div_item = tsne_plot(outdot,outdiv,outhtml)
    # high rep gene in each cluster
    top_gene_cluster(matfile,tsnefile,outtbl,args.mincell,args.topgene)
    tbltxt = pd.read_csv(outtbl,header=0,sep="\t")
    tbldiv = tsne_table(tbltxt)
    #with open(outdiv,'a') as w:
    #    w.write('\n'+tbldiv)

    # report
    report_prepare(div_item, tbldiv, step)
    step.clean_up()



def get_opts_replace_tsne(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument('--tsne', help='tsne file', required=True)
        parser.add_argument('--mat', help='matrix rep file', required=True)
        parser.add_argument('--rep', help='cell rep file', required=True)
        parser.add_argument('--mincell', type=int, default=5, help='turn-over in at least cells, default 5')
        parser.add_argument('--topgene', type=int, default=10, help='top N genes,default 10')


