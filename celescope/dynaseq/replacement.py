#!/bin/env python
# coding=utf8

import os
import sys
import subprocess
import pandas as pd
import pysam
from celescope.tools.step import Step, s_common
import celescope.tools.utils as utils

toolsdir = os.path.dirname(__file__)


class Replacement(Step):
    """
    Features
    - Computes the replacement rates in each cell and gene.
    - Boxplots for rates distribution.

    Output
    - `{sample}.TC_matrix.rds` New and old info for each barcode/gene/umi.
    - `{sample}.new_matrix.tsv.gz` New RNA matrix.
    - `{sample}.old_matrix.tsv.gz` Old RNA matrix.
    - `{sample}.fraction_of_newRNA_per_cell.txt` Fraction of new RNA of each cell.
    - `{sample}.fraction_of_newRNA_per_gene.txt` Fraction of new RNA of each gene.
    - `{sample}.fraction_of_newRNA_matrix.txt` Fraction of new RNA of each cell and gene.
    """

    def __init__(self, args):
        Step.__init__(self, args)

        # input files
        self.outdir = args.outdir
        self.sample = args.sample
        self.bam_file = args.bam
        self.snp_file = args.bg
        self.bg_cov = args.bg_cov
        self.cell_keep = args.cell_keep
        # output files
        self.outread = os.path.join(self.outdir, self.sample+'.corrected_gene_cell_UMI_read.txt')
        self.outrds = os.path.join(self.outdir, self.sample+'.TC_matrix.rds')
        self.outpre = os.path.join(self.outdir, self.sample)

    @utils.add_log
    def run(self):
        # get backgroud snp
        bg = self.background_snp(self.snp_file, self.bg_cov)
        # get reads with TC
        self.extract_dem(self.bam_file, self.outread, bg)
        # run_R
        self.generate_TC_matrix(self.outread, self.outrds, self.cell_keep)

        # split to New and Old Matrix
        totMat = self.outrds+'.tsv'
        new_mat = self.outpre+'.new_matrix.tsv'
        old_mat = self.outpre+'.old_matrix.tsv'
        con_mat = self.outpre+'.NvsO_matrix.tsv'
        self.split_matrix(totMat, self.outpre)

        # replacement stat
        self.replacment_stat(con_mat, self.outpre)
        # plot
        div_item = self.replacment_plot(self.outpre)

        # report
        self.report_prepare(div_item)
        self._clean_up()

        # clean
        cmd = ['rm', self.outread]
        self.run_cmd(cmd)
        cmd = ['rm', self.outrds+'.tsv']
        self.run_cmd(cmd)
        cmd = ['rm', con_mat]
        self.run_cmd(cmd)
        cmd = ['gzip', new_mat]
        self.run_cmd(cmd)
        cmd = ['gzip', old_mat]
        self.run_cmd(cmd)

    def run_cmd(self, cmd):
        subprocess.call(' '.join(cmd), shell=True)

    @utils.add_log
    def extract_dem(self, bam, outfile, bg):
        bamfile = pysam.AlignmentFile(bam, 'rb')
        countdict = {}
        for read in bamfile.fetch():
            try:
                chro = read.reference_name
                cb = read.get_tag('CB')
                ub = read.get_tag('UB')
                if not read.has_tag('GN'):
                    continue
                gene = read.get_tag('GN')

                if read.get_tag('ST') == '+':
                    stag = read.get_tag('TL')
                else:
                    stag = read.get_tag('AL')
                if len(stag) == 1 and stag[0] == 0:
                    gene += '--T'
                else:
                    fcount = 0
                    for si in range(0, len(stag)):
                        pos = chro + '_' + str(stag[si])
                        if pos in bg:
                            fcount += 1
                    if fcount == len(stag):
                        gene += '--T'
                    else:
                        gene += '--C'

                readinfo = '\t'.join([gene, cb, ub])
                if readinfo not in countdict:
                    countdict[readinfo] = 1
                else:
                    countdict[readinfo] += 1

            except (ValueError, KeyError):
                continue
        bamfile.close()

        out1 = open(outfile, 'w')
        for rid in countdict:
            out1.write(rid+'\t'+str(countdict[rid])+'\n')
        out1.close()

    @utils.add_log
    def background_snp(self, bgfile, cov=1):
        outdict = {}
        if bgfile.endswith('.csv'):
            with open(bgfile) as f:
                f.readline()
                for i in f:
                    ii = i.strip().split(',')
                    if int(ii[2]) < cov:
                        continue
                    chr_pos = ii[4]+'_'+ii[1]
                    outdict[chr_pos] = 1
        elif bgfile.endswith('.vcf'):
            bcf_in = pysam.VariantFile(bgfile)
            for rec in bcf_in.fetch():
                try:
                    chrom, pos = rec.chrom, rec.pos
                    chr_pos = chrom+'_'+str(pos)
                    outdict[chr_pos] = 1
                except (ValueError, KeyError):
                    continue
            bcf_in.close()

        else:
            try:
                sys.exit(1)
            except SystemExit:
                print('Background snp file format cannot be recognized! Only csv or vcf format.')
            finally:
                print('Background snp file format cannot be recognized! Only csv or vcf format.')
        return outdict

    @utils.add_log
    def generate_TC_matrix(self, read, outrds, cell=100000):
        app = toolsdir + "/Generate_T_C_matrix.R"
        cmd = (
            f'Rscript {app} {read} {cell} {outrds}'
        )
        os.system(cmd)

    @utils.add_log
    def split_matrix(self, mat, outpre):
        outnew = open(outpre+'.new_matrix.tsv', 'w')
        outold = open(outpre+'.old_matrix.tsv', 'w')
        con_mat = open(outpre+'.NvsO_matrix.tsv', 'w')
        infile = open(mat, 'r')

        tmph = infile.readline().strip().split()
        fill_na = ['0'] * len(tmph)
        tmph.insert(0, '')
        outnew.write('\t'.join(tmph)+'\n')
        outold.write('\t'.join(tmph)+'\n')
        con_mat.write('\t'.join(tmph)+'\n')

        genes = {}
        for i in infile:
            ii = i.strip().split()
            gt = ii[0].split('--')
            ii[0] = gt[0]
            if gt[0] not in genes:
                genes[gt[0]] = [0, [], []]
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
            if genes[gi][0] == 3:
                for ci in range(len(genes[gi][1])):
                    con_mat.write('\t'+genes[gi][1][ci]+':'+genes[gi][2][ci])
                con_mat.write('\n')
            elif genes[gi][0] == 1:
                outold.write(gi+'\t'+'\t'.join(fill_na)+'\n')
                for ci in range(len(genes[gi][1])):
                    con_mat.write('\t'+genes[gi][1][ci]+':'+'0')
                con_mat.write('\n')
            elif genes[gi][0] == 2:
                outnew.write(gi+'\t'+'\t'.join(fill_na)+'\n')
                for ci in range(len(genes[gi][2])):
                    con_mat.write('\t'+'0'+':'+genes[gi][2][ci])
                con_mat.write('\n')

        outnew.close()
        outold.close()
        con_mat.close()
        infile.close()

    @utils.add_log
    def replacment_stat(self, inmat, outpre, mincell=10, mingene=10):

        outcell = open(outpre+'.fraction_of_newRNA_per_cell.txt', 'w')
        outgene = open(outpre+'.fraction_of_newRNA_per_gene.txt', 'w')
        outmat = open(outpre+'.fraction_of_newRNA_matrix.txt', 'w')

        cells = {}
        genes = {}
        mats = {}
        with open(inmat) as f:
            hh = f.readline().strip().split()
            hh.insert(0, '')
            outmat.write('\t'.join(hh)+'\n')
            for h in hh[1:]:
                cells[h] = [[], []]
            for i in f:
                ii = i.strip().split()
                mats[ii[0]] = []
                genes[ii[0]] = [[], []]
                for xi in range(1, len(ii)):
                    xx = [int(x) for x in ii[xi].split(':')]
                    if sum(xx) == 0:
                        tmpf = 'NA'
                    else:
                        tmpf = float(xx[0])/(xx[0]+xx[1])
                    mats[ii[0]].append(str(tmpf))
                    if sum(xx) < 2:
                        continue
                    cells[hh[xi]][0].append(float(xx[0]))
                    cells[hh[xi]][1].append(xx[1])
                    genes[ii[0]][0].append(float(xx[0]))
                    genes[ii[0]][1].append(xx[1])

        for ci in cells:
            if len(cells[ci][0]) < mincell:
                continue
            cfloat = sum(cells[ci][0])/(sum(cells[ci][0])+sum(cells[ci][1]))
            outcell.write(ci+'\t'+str(cfloat)+'\n')

        for gi in genes:
            if len(genes[gi][0]) < mingene:
                continue
            gfloat = sum(genes[gi][0])/(sum(genes[gi][0])+sum(genes[gi][1]))
            outgene.write(gi+'\t'+str(gfloat)+'\n')

        for mi in mats:
            outmat.write(mi+'\t'+'\t'.join(mats[mi])+'\n')

        outcell.close()
        outgene.close()
        outmat.close()

    @utils.add_log
    def replacment_plot(self, sample):
        import plotly
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        outpre = os.path.basename(sample)
        df1 = pd.read_table(sample+'.fraction_of_newRNA_per_gene.txt', header=None)
        df1.columns = ['gene', 'per']
        df = pd.read_table(sample+'.fraction_of_newRNA_per_cell.txt', header=None)
        df.columns = ['gene', 'per']

        fig = make_subplots(rows=1, cols=2)
        fig.add_trace(
            go.Violin(y=df1['per'], box_visible=True, line_color='black',
                      meanline_visible=True, fillcolor='#1f77b4', opacity=0.6, x0=outpre),
            row=1, col=1
        )
        fig.add_trace(
            go.Violin(y=df['per'], box_visible=True, line_color='black',
                      meanline_visible=True, fillcolor='#ff7f0e', opacity=0.6, x0=outpre),
            row=1, col=2
        )

        fig.update_layout(yaxis_zeroline=True,  showlegend=False)
        fig.update_layout(plot_bgcolor='#FFFFFF')
        fig.update_xaxes(showgrid=False, linecolor='black', showline=True, ticks=None)
        fig.update_yaxes(showgrid=False, linecolor='black', showline=True, ticks='outside',
                         title_text="Fraction new RNA (per gene)", row=1, col=1, rangemode="tozero")
        fig.update_yaxes(showgrid=False, linecolor='black', showline=True, ticks='outside',
                         title_text="Fraction new RNA (per cell)", row=1, col=2, rangemode="tozero")

        div = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

        return div

    def report_prepare(self, outdiv):
        self.add_data(replacement=outdiv)


@utils.add_log
def replacement(args):

    with Replacement(args) as runner:
        runner.run()


def get_opts_replacement(parser, sub_program):
    parser.add_argument('--bg_cov', type=int, default=1,
                        help='background snp depth filter, lower than bg_cov will be discarded. Only valid in csv format')
    if sub_program:
        parser.add_argument('--bam', help='bam file from conversion step', required=True)
        parser.add_argument('--bg', help='background snp file, csv or vcf format', required=True)
        parser.add_argument('--cell_keep', type=int, default=100000, help='filter cell')
        parser.add_argument('--min_cell', type=int, default=10, help='a gene expressed in at least cells, default 10')
        parser.add_argument('--min_gene', type=int, default=10, help='at least gene num in a cell, default 10')
        parser = s_common(parser)
    return parser
