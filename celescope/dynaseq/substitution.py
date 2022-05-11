#!/bin/env python
# coding=utf8

import os
import pysam
import re
import pandas as pd
import plotly
import plotly.graph_objects as go
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT

class Substitution(Step):
    """
    ## Features
    - Computes the overall conversion rates in reads and plots a barplot.

    ## Output
    - `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.
    """

    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        # input files
        self.bam_file = args.bam

        # output files
        self.outstat = os.path.join(self.outdir, self.sample+'.substitution.txt')

    @utils.add_log
    def run(self):
        # overall rate
        for_base, rev_base, is_forward, is_reverse = self.get_sub_tag(self.bam_file)
        self.sub_stat(for_base, rev_base, is_forward, is_reverse, self.outstat)
        div_item = self.sub_plot(self.outstat)

        self.report_prepare(div_item)

    @utils.add_log
    def get_sub_tag(self, bam):
        bamfile = pysam.AlignmentFile(bam, 'rb')
        is_reverse = {'cA': 0, 'gA': 0, 'tA': 0, 'aC': 0, 'gC': 0,
                      'tC': 0, 'aG': 0, 'cG': 0, 'tG': 0, 'aT': 0, 'cT': 0, 'gT': 0}
        is_forward = {'cA': 0, 'gA': 0, 'tA': 0, 'aC': 0, 'gC': 0,
                      'tC': 0, 'aG': 0, 'cG': 0, 'tG': 0, 'aT': 0, 'cT': 0, 'gT': 0}
        for_base = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        rev_base = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        snp_tags = ['', 'cA', 'gA', 'tA', 'aC', 'gC', 'tC', 'aG', 'cG', 'tG', 'aT', 'cT', 'gT']
        ref_tags = ['', 'a', 'c', 'g', 't']
        for read in bamfile.fetch():
            try:
                snpmatch = re.match(
                    r'cA(\d+);gA(\d+);tA(\d+);aC(\d+);gC(\d+);tC(\d+);aG(\d+);cG(\d+);tG(\d+);aT(\d+);cT(\d+);gT(\d+);', read.get_tag('SC'), re.M)
                totmatch = re.match(r'a(\d+);c(\d+);g(\d+);t(\d+)', read.get_tag('TC'), re.M)
                if snpmatch and totmatch:
                    if read.is_reverse:
                        for j in range(1, len(ref_tags)):
                            rev_base[ref_tags[j]] += int(totmatch.group(j))
                        for i in range(1, len(snp_tags)):
                            is_reverse[snp_tags[i]] += int(snpmatch.group(i))
                    else:
                        for j in range(1, len(ref_tags)):
                            for_base[ref_tags[j]] += int(totmatch.group(j))
                        for i in range(1, len(snp_tags)):
                            is_forward[snp_tags[i]] += int(snpmatch.group(i))
            except (ValueError, KeyError):
                continue
        bamfile.close()

        return for_base, rev_base, is_forward, is_reverse

    @utils.add_log
    def sub_stat(self, for_base, rev_base, is_forward, is_reverse, outfile):
        convertdict = {'a': ['aC', 'aG', 'aT'],
                       'c': ['cA', 'cG', 'cT'],
                       'g': ['gA', 'gC', 'gT'],
                       't': ['tA', 'tC', 'tG']}
        subdict = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                   'aC': 'tG', 'aG': 'tC', 'aT': 'tA',
                   'cA': 'gT', 'cG': 'gC', 'cT': 'gA',
                   'gA': 'cT', 'gC': 'cG', 'gT': 'cA',
                   'tA': 'aT', 'tC': 'aG', 'tG': 'aC'}
        outdict = {'aC': 'A_to_C', 'aG': 'A_to_G', 'aT': 'A_to_T',
                   'cA': 'C_to_A', 'cG': 'C_to_G', 'cT': 'C_to_T',
                   'gA': 'G_to_A', 'gC': 'G_to_C', 'gT': 'G_to_T',
                   'tA': 'T_to_A', 'tC': 'T_to_C', 'tG': 'T_to_G'}
        outw = open(outfile, 'w')
        for x in ['a', 'c', 'g', 't']:
            fbase = for_base[x]
            rbase = rev_base[subdict[x]]
            for y in convertdict[x]:
                fcov = is_forward[y]*100 / float(fbase)
                rcov = is_reverse[subdict[y]]*100 / float(rbase)
                outw.write(outdict[y]+'\t'+"%.3f" % fcov+'\t'+"%.3f" % rcov+'\n')
        outw.close()

    @utils.add_log
    def sub_plot(self, txt):
        df = pd.read_table(txt, header=None)
        df.columns = ['sample', '+', '-']

        fig = go.Figure()
        # 设置颜色：
        import plotly.express as px
        num4colors = 0
        num4rainbow = 0
        colors_list = []
        while num4colors < 100:
            if num4rainbow == 9:
                num4rainbow = 0
            colors_list.append(px.colors.qualitative.Plotly[num4rainbow])
            num4colors += 1
            num4rainbow += 1

        num4sample = 0
        colors4sample = {}
        num4x = 0

        for sample in df['sample'].unique():
            legend_show = True
            colors4sample[sample] = colors_list[num4sample]
            num4sample += 1
            flag_x = 'x' + str(num4x+1)
            df_plot = df[df['sample'] == sample]
            num4x += 1

            fig.add_trace(go.Bar(name=sample+'+',
                                 x=df_plot['sample'],
                                 y=df_plot['+'],
                                 legendgroup=sample,
                                 marker_color=colors4sample[sample],
                                 marker_line_color='#FFFFFF',
                                 showlegend=legend_show,
                                 xaxis=flag_x)
                          )
            fig.add_trace(go.Bar(name=sample+'-',
                                 x=df_plot['sample'],
                                 y=df_plot['-'],
                                 legendgroup=sample,
                                 showlegend=legend_show,
                                 marker_color=colors4sample[sample],
                                 marker_line_color='#FFFFFF',
                                 opacity=0.3,
                                 xaxis=flag_x)
                          )

        fig.update_layout(barmode='stack')

        per = 1/(num4x+1)
        gap4bar = per/len(df['sample'].unique())
        num4x = 0
        for typeB in df['sample'].unique():
            if num4x == 0:
                flag_x = 'xaxis'
            else:
                flag_x = 'xaxis' + str(num4x+1)
            anchor_x = 'x'+str(num4x+1)
            num4x += 1
            fig['layout'][flag_x] = dict(domain=[per*num4x, per*(num4x+1)-gap4bar], anchor=anchor_x, title=typeB)

        fig.update_layout(plot_bgcolor='#FFFFFF')
        fig.update_xaxes(showgrid=False, linecolor='black', showline=True, ticks='outside', showticklabels=False)
        fig.update_yaxes(showgrid=False, linecolor='black', showline=True, ticks='outside')
        width_num = 400 * (len(df['sample'].unique()) * len(df['sample'].unique())) / (5*12)  # 控制柱形图的宽度
        fig.update_layout(height=500, width=width_num)
        fig.update_layout(legend=dict(orientation="h"))
        fig.update_layout(legend=dict(
            yanchor="top",
            y=1.3,
            xanchor="left",
            x=0.05,
            valign="top",
        ))

        fig.update_layout(
            yaxis_title="Rates of nucleotide substitution (%)",
        )
        fig.update_xaxes(
            tickangle=-80,
            title_font={"size": 15},
            title_standoff=25
        )

        div = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

        return div

    def report_prepare(self, outdiv):
        self.add_data(substitution=outdiv)


@utils.add_log
def substitution(args):

    with Substitution(args,display_title='Substitution') as runner:
        runner.run()


def get_opts_substitution(parser, sub_program):
    if sub_program:
        parser.add_argument('--bam', help=HELP_DICT['bam_for_replacement'], required=True)
        parser = s_common(parser)
    return parser
