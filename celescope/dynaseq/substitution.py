#!/bin/env python
# coding=utf8

import os
import pysam
import re
import pandas as pd
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.tools.plotly_plot import Substitution_plot



class Substitution(Step):
    """
    ## Features
    - Computes the overall conversion rates in reads and plots a barplot.

    ## Output
    - `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.
    """

    def __init__(self, args):
        Step.__init__(self, args)

        # input files
        self.sample = args.sample
        self.bam_file = args.bam
        self.outdir = args.outdir

        # output files
        self.outstat = os.path.join(self.outdir, self.sample+'.substitution.txt')

    @utils.add_log
    def run(self):
        # overall rate
        for_base, rev_base, is_forward, is_reverse = self.get_sub_tag(self.bam_file)
        self.sub_stat(for_base, rev_base, is_forward, is_reverse, self.outstat)
        #self.add_substitution_metrics()
        self.substitution_plot()
        self.add_help()

    @utils.add_log
    def get_sub_tag(self, bam):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, 'rb', require_index=False)
        pysam.set_verbosity(save)
        is_reverse = {'cA': 0, 'gA': 0, 'tA': 0, 'aC': 0, 'gC': 0,
                      'tC': 0, 'aG': 0, 'cG': 0, 'tG': 0, 'aT': 0, 'cT': 0, 'gT': 0}
        is_forward = {'cA': 0, 'gA': 0, 'tA': 0, 'aC': 0, 'gC': 0,
                      'tC': 0, 'aG': 0, 'cG': 0, 'tG': 0, 'aT': 0, 'cT': 0, 'gT': 0}
        for_base = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        rev_base = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        snp_tags = ['', 'cA', 'gA', 'tA', 'aC', 'gC', 'tC', 'aG', 'cG', 'tG', 'aT', 'cT', 'gT']
        ref_tags = ['', 'a', 'c', 'g', 't']
        
        for read in bamfile.fetch(until_eof=True):
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
    def substitution_plot(self):
        df = pd.read_table(self.outstat, header=None)
        df.columns = ['sample', '+', '-']
        subplot = Substitution_plot(df_bar=df).get_plotly_div()
        self.add_data(substitution_box=subplot)

    @utils.add_log
    def add_help(self):
        self.add_help_content(
            name='Substitution:',
            content='Nucleotide substitution rate per conversion type'
        )


@utils.add_log
def substitution(args):

    with Substitution(args) as runner:
        runner.run()


def get_opts_substitution(parser, sub_program):
    if sub_program:
        parser.add_argument('--bam', help='bam file from conversion step', required=True)
        parser = s_common(parser)
    return parser
