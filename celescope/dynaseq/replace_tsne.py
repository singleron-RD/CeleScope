#!/bin/env python
# coding=utf8

import os
import pandas as pd
import plotly
import plotly.graph_objects as go
from celescope.tools.step import Step, s_common
from celescope.tools import utils


class Replace_tsne(Step):
    """
    ## Features
    - Replace rate in each cluster
    - Top replace genes in each cluster

    ## Output
    - `{sample}.rep_in_tsne.txt` Replace rate in each cluster.
    - `{sample}.rep_in_tsne_top10` Top 10 replace genes in each cluster.
    """

    def __init__(self, args):
        Step.__init__(self, args)

        # input files
        self.sample = args.sample
        self.tsnefile = args.tsne
        self.matfile = args.mat
        self.repfile = args.rep
        self.mincell = args.mincell
        self.topgene = args.topgene
        # output files
        self.outdot = os.path.join(self.outdir, self.sample+'.rep_in_tsne.txt')
        self.outtbl = os.path.join(self.outdir, self.sample+'.rep_in_tsne_top10.txt')

    @utils.add_log
    def run(self):
        # rep in cells in cluster
        self.dot_tsne(self.repfile, self.tsnefile, self.outdot)
        div_item = self.tsne_plot(self.outdot)
        # high rep gene in each cluster
        self.top_gene_cluster(self.matfile, self.tsnefile, self.outtbl, self.mincell, self.topgene)
        tbltxt = pd.read_csv(self.outtbl, header=0, sep="\t")
        tbldiv = self.tsne_table(tbltxt)

        # report
        self.report_prepare(div_item, tbldiv)
        self._clean_up()

    @utils.add_log
    def dot_tsne(self, repfile, tsnefile, outfile):
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

    @utils.add_log
    def tsne_plot(self, txt):
        df = pd.read_table(txt)
        df.sort_values(by="ratio")
        newtitle = "t-SNE plot Colored by RNA Turn-over rate"

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=df['tSNE_1'], y=df['tSNE_2'], mode='markers',
                                 marker_opacity=0.9, marker_size=4, marker_color=df['ratio'],
                                 marker_colorscale="PuBu", marker_showscale=True,
                                 ))
        fig.update_layout(height=600, width=600, title_text=newtitle)
        fig.update_layout(plot_bgcolor='#FFFFFF')
        fig.update_xaxes(showgrid=False, linecolor='black', showline=True, ticks='outside', title_text='t-SNE1')
        fig.update_yaxes(showgrid=False, linecolor='black', showline=True, ticks='outside', title_text='t-SNE2')

        div = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

        return div

    def tsne_table(self, txt):
        marker_gene_table = txt.to_html(
            escape=False,
            index=False,
            table_id='replacement_table_cluster',
            justify="center")

        return marker_gene_table

    def file_stat(self, infile, clu):
        clus = list(set(clu.values()))
        cluster = {}
        for c in clus:
            cluster[c] = {}
        fn = open(infile, "r")
        fnh = fn.readline().strip().split()
        fnh.insert(0, '')
        for i in fn:
            ii = i.strip().split()
            for j in range(1, len(ii)):
                if ii[j] == 'NA':
                    continue
                if fnh[j] not in clu:
                    continue
                if ii[0] not in cluster[clu[fnh[j]]]:
                    cluster[clu[fnh[j]]][ii[0]] = []
                cluster[clu[fnh[j]]][ii[0]].append(float(ii[j]))
        fn.close()
        return cluster

    def tsne_file(self, infile):
        clu = {}
        with open(infile) as f:
            f.readline()
            for i in f:
                ii = i.strip().split()
                clu[ii[0]] = ii[3]
        return clu

    @utils.add_log
    def top_gene_cluster(self, matrix, tsnefile, outfile, mincell=5, topgene=10):
        tsne = self.tsne_file(tsnefile)
        cluster = self.file_stat(matrix, tsne)

        w = open(outfile, 'w')
        w.write("cluster\tgene\tTurn-over_rate\tcells\n")
        for c in cluster:
            tmp = {}
            for g in cluster[c]:
                gt = sum(cluster[c][g]) / len(cluster[c][g])
                tmp[g] = gt
            sorttmp = sorted(tmp.items(), key=lambda item: item[1], reverse=True)
            tmpn = 0
            for x in sorttmp:
                if len(cluster[c][x[0]]) < mincell:
                    continue
                tmpn += 1
                if tmpn > topgene:
                    break
                w.write('cluster'+c+'\t'+x[0]+'\t'+str(x[1])+'\t'+str(len(cluster[c][x[0]]))+'\n')
        w.close()

    def report_prepare(self, outdiv, outable):
        self.add_data(replace_tsne=outdiv)
        self.add_data(replace_tsne_table=outable)


@utils.add_log
def replace_tsne(args):

    with Replace_tsne(args) as runner:
        runner.run()


def get_opts_replace_tsne(parser, sub_program):
    if sub_program:
        parser.add_argument('--tsne', help='tsne file from analysis step', required=True)
        parser.add_argument('--mat', help='matrix replacement file, from replacement step', required=True)
        parser.add_argument('--rep', help='cell replacement file, from replacement step', required=True)
        parser.add_argument('--mincell', type=int, default=5, help='turn-over in at least cells, default 5')
        parser.add_argument('--topgene', type=int, default=10, help='show top N genes,default 10')
        parser = s_common(parser)
    return parser
