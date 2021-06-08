#!/bin/env python
# coding=utf8

import os
import sys
import argparse
import json
import logging
import subprocess
import re
import numpy as np
import pandas as pd
import glob
import pysam
#from celescope.tools.report import reporter
from celescope.tools.utils import *
from celescope.tools.Step import Step, s_common

toolsdir = os.path.dirname(__file__)


@add_log
def extract_dem(bam,outfile,bg):
    bamfile = pysam.AlignmentFile(bam, 'rb')
    countdict = {}
    for read in bamfile.fetch():
        try:
            chr = read.reference_name
            cb = read.get_tag('CB')
            ub = read.get_tag('UB')
            if not read.has_tag('GN'): continue
            gene = read.get_tag('GN')

            if read.get_tag('ST') == '+':
                stag = read.get_tag('TL')
            else:
                stag = read.get_tag('AL')
            if len(stag)==1 and stag[0]==0:
                gene += '--T'
            else:
                fcount = 0
                for si in range(0,len(stag)):
                    pos = chr + '_' + str(stag[si])
                    if pos in bg:
                        fcount += 1
                if fcount == len(stag):
                    gene += '--T'
                else:
                    gene += '--C'
            
            readinfo = '\t'.join([gene,cb,ub])
            if readinfo not in countdict:
                countdict[readinfo] = 1
            else:
                countdict[readinfo] += 1

        except (ValueError,KeyError):
            continue
    bamfile.close()

    print('Write gene_cell_UMI_read to {}'.format(outfile))
    out1 = open(outfile,'w')
    for rid in countdict:
        out1.write(rid+'\t'+str(countdict[rid])+'\n')
    out1.close()

@add_log
def background_snp(bgfile,cov=1):
    outdict = {}
    if bgfile.endswith('.csv'):
        with open(bgfile) as f:
            f.readline()
            for i in f:
                ii = i.strip().split(',')
                if int(ii[2])<cov: continue
                chr_pos = ii[1]+'_'+ii[5]
                outdict[chr_pos] = 1
    elif bgfile.endswith('.vcf'):
        print('TBC')
    else:
        try:
            sys.exit(1)
        except SystemExit:
            print('Background snp file format cannot be recognized! Only csv and vcf format.')
        finally:
            print('Background snp file format cannot be recognized! Only csv and vcf format.')
    return outdict





    """
    return html code
    """
    marker_df = marker_df.loc[:, ["cluster", "gene",
                                  "avg_logFC", "pct.1", "pct.2", "p_val_adj"]]
    marker_df["cluster"] = marker_df["cluster"].apply(lambda x: f"cluster {x}")
    marker_gene_table = marker_df.to_html(
        escape=False,
        index=False,
        table_id="marker_gene_table",
        justify="center")
    return marker_gene_table


@add_log
def generate_TC_matrix(read, outrds, cell=100000):
    app = toolsdir + "/Generate_T_C_matrix.R"
    cmd = (
        f'Rscript {app} {read} {cell} {outrds}'
    )
    os.system(cmd)

@add_log
def split_matrix(mat,outpre):
    outnew = open(outpre+'.new_matrix.tsv', 'w')
    outold = open(outpre+'.old_matrix.tsv', 'w')
    con_mat = open(outpre+'.NvsO_matrix.tsv', 'w')
    infile = open(mat, 'r')

    tmph = infile.readline().strip().split()
    fill_na = ['0'] * len(tmph)
    tmph.insert( 0, 'geneID')
    outnew.write('\t'.join(tmph)+'\n')
    outold.write('\t'.join(tmph)+'\n')
    con_mat.write('\t'.join(tmph)+'\n')

    genes = {}
    for i in infile:
        ii = i.strip().split()
        gt = ii[0].split('--')
        ii[0] = gt[0]
        if gt[0] not in genes:
            genes[gt[0]] = [0,[],[]]
        if gt[1] == 'C':
            genes[gt[0]][0] += 1
            genes[gt[0]][1] = ii[1:]
            outnew.write('\t'.join(ii)+'\n')
        elif gt[1] == 'T':
            genes[gt[0]][0] += 2
            genes[gt[0]][2] = ii[1:]
            outold.write('\t'.join(ii)+'\n')
    
    for gi in genes:
        con_mat.write(gi)
        if genes[gi][0]==3:
            for ci in range(len(genes[gi][1])):
                con_mat.write('\t'+genes[gi][1][ci]+':'+genes[gi][2][ci])
            con_mat.write('\n')
        elif genes[gi][0]==1:
            outold.write(gi+'\t'+'\t'.join(fill_na)+'\n')
            for ci in range(len(genes[gi][1])):
                con_mat.write('\t'+genes[gi][1][ci]+':'+'0')
            con_mat.write('\n')
        elif genes[gi][0]==2:
            outnew.write(gi+'\t'+'\t'.join(fill_na)+'\n')
            for ci in range(len(genes[gi][2])):
                con_mat.write('\t'+'0'+':'+genes[gi][2][ci])
            con_mat.write('\n')

    outnew.close()
    outold.close()
    con_mat.close()
    infile.close()



@add_log
def replacment_stat(inmat,outpre,mincell=10,mingene=10):

    outcell = open(outpre+'.fraction_of_newRNA_per_cell.txt', 'w')
    outgene = open(outpre+'.fraction_of_newRNA_per_gene.txt', 'w')
    outmat = open(outpre+'.fraction_of_newRNA_matrix.txt', 'w')

    cells = {}
    genes = {}
    mats = {}
    with open(inmat) as f:
        hh = f.readline().strip().split()
        outmat.write('\t'.join(hh)+'\n')
        for h in hh[1:]:
            cells[h] = [[],[]]
        for i in f:
            ii = i.strip().split()
            mats[ii[0]] = []
            genes[ii[0]] = [[],[]]
            for xi in range(1,len(ii)):
                xx = [int(x) for x in ii[xi].split(':')]                 
                if sum(xx) == 0:
                    tmpf = 'NA'                   
                else:
                    tmpf = float(xx[0])/(xx[0]+xx[1])
                mats[ii[0]].append(str(tmpf))
                if sum(xx)<2:
                    continue                
                cells[hh[xi]][0].append(float(xx[0]))
                cells[hh[xi]][1].append(xx[1])
                genes[ii[0]][0].append(float(xx[0]))
                genes[ii[0]][1].append(xx[1])

    for ci in cells:
        if len(cells[ci][0])<mincell: continue
        cfloat = sum(cells[ci][0])/(sum(cells[ci][0])+sum(cells[ci][1]))
        outcell.write(ci+'\t'+str(cfloat)+'\n')
    
    for gi in genes:
        if len(genes[gi][0])<mingene: continue
        gfloat = sum(genes[gi][0])/(sum(genes[gi][0])+sum(genes[gi][1]))
        outgene.write(gi+'\t'+str(gfloat)+'\n')

    for mi in mats:
        outmat.write(mi+'\t'+'\t'.join(mats[mi])+'\n')

    outcell.close()
    outgene.close()
    outmat.close()


@add_log
def replacment_plot(sample):
    import pandas as pd
    import plotly
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    outpre = os.path.basename(sample)
    outdir = os.path.dirname(sample)
    df1 = pd.read_table(sample+'.fraction_of_newRNA_per_gene.txt', header=None)
    df1.columns = ['gene', 'per']
    df = pd.read_table(sample+'.fraction_of_newRNA_per_cell.txt', header=None)
    df.columns = ['gene', 'per']

    fig = make_subplots(rows=1, cols=2)
    fig.add_trace(
        go.Violin(y=df1['per'], box_visible=True, line_color='black',
              meanline_visible=True, fillcolor='#1f77b4', opacity=0.6, x0=outpre) ,
        row=1, col=1
    )
    fig.add_trace(
        go.Violin(y=df['per'], box_visible=True, line_color='black',
              meanline_visible=True, fillcolor='#ff7f0e', opacity=0.6, x0=outpre) ,
        row=1, col=2
    )

    fig.update_layout(yaxis_zeroline=True,  showlegend=False)
    fig.update_layout(plot_bgcolor = '#FFFFFF')
    fig.update_xaxes(showgrid=False, linecolor='black', showline=True, ticks=None)
    fig.update_yaxes(showgrid=False, linecolor='black', showline=True, ticks='outside',title_text="Fraction new RNA (per gene)",row=1, col=1, rangemode="tozero")
    fig.update_yaxes(showgrid=False, linecolor='black', showline=True, ticks='outside',title_text="Fraction new RNA (per cell)",row=1, col=2, rangemode="tozero")

    div = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')
    '''
    with open(sample+'.replacement.div', 'w') as fno:
        fno.write("<h3>RNA Turn-over rate</h3>\n")
        fno.write(div + '\n')
    fig.write_html(sample+".replacement.html")
    '''
    return div

def report_prepare(outdiv, step_obj):
    step_obj.add_data_item(replacement=outdiv)


@add_log
def replacement(args):

    step_name = "replacement"
    step = Step(args, step_name)

    # check dir
    outdir = args.outdir
    sample = args.sample
    assay = args.assay 

    bam_file = args.bam
    snp_file = args.bg

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))
    outread = os.path.join(outdir, sample+'.corrected_gene_cell_UMI_read.txt')
    outrds = os.path.join(outdir, sample+'.TC_matrix.rds')
    outpre = os.path.join(outdir,sample)
    
    # get backgroud snp
    bg = background_snp(snp_file,args.bg_cov)
    print('Parse backgroup pos: {}'.format(str(len(bg))))

    # get reads with TC
    extract_dem(bam_file,outread,bg)

    # run_R
    generate_TC_matrix(outread, outrds, args.cell_keep)
    print('Write matrix to : {}'.format(outrds))

    # split to New and Old Matrix
    totMat = outrds+'.tsv'
    new_mat = outpre+'.new_matrix.tsv'
    old_mat = outpre+'.old_matrix.tsv'
    con_mat = outpre+'.NvsO_matrix.tsv'
    split_matrix(totMat,outpre)
    print('Write New RNA matrix to : {}'.format(new_mat))
    print('Write Old RNA matrix to : {}'.format(old_mat))
    
    # replacement stat
    replacment_stat(con_mat,outpre)

    # plot
    div_item = replacment_plot(outpre)

    # report
    report_prepare(div_item, step)
    step.clean_up()

    # clean
    def run_cmd(cmd):
        subprocess.call(' '.join(cmd),shell=True)
    print('clean tmp files ...')
    cmd=['rm', outread]
    run_cmd(cmd)
    cmd=['rm', outrds+'.tsv']
    run_cmd(cmd)
    cmd=['rm', con_mat]
    run_cmd(cmd)
    cmd=['gzip', new_mat]
    run_cmd(cmd)
    cmd=['gzip', old_mat]
    run_cmd(cmd)


def get_opts_replacement(parser, sub_program):
    parser.add_argument('--bg_cov', type=int, default=1, help='background snp depth filter, lower than bg_cov will be discarded ')
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--bam', help='bam file', required=True)
        parser.add_argument('--bg', help='background snp file', required=True)
        parser.add_argument('--cell_keep', type=int, default=100000, help='filter cell')
        parser.add_argument('--min_cell', type=int, default=10, help='a gene expressed in at least cells, default 10')
        parser.add_argument('--min_gene', type=int, default=10, help='at least gene num in a cell, default 10')


