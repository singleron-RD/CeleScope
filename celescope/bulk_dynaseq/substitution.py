#!/bin/env python
# coding=utf8

import os
import pysam
import re
import pandas as pd
from multiprocessing import Pool
from celescope.tools.plotly_plot import Substitution_plot_bulk
from celescope.tools.step import Step, s_common
from celescope.tools import utils

sorted_header = [
    "A_to_C",
    "A_to_G",
    "A_to_T",
    "C_to_A",
    "C_to_G",
    "C_to_T",
    "G_to_A",
    "G_to_C",
    "G_to_T",
    "T_to_A",
    "T_to_C",
    "T_to_G",
]
tc_header = [
    "T_to_C",
    "A_to_C",
    "A_to_G",
    "A_to_T",
    "C_to_A",
    "C_to_G",
    "C_to_T",
    "G_to_A",
    "G_to_C",
    "G_to_T",
    "T_to_A",
    "T_to_G",
]


class Substitution(Step):
    """
    ## Features
    - Computes the overall conversion rates in reads per well.
    - Per-well outputs (no merging).

    ## Output
    - `{well}/{well}.substitution.txt` Per-well overall conversion rates.
    - `{well}/{well}.TC_substitution.txt` Per-well TC substitution rates.
    """

    def __init__(self, args):
        Step.__init__(self, args)

        # set
        self.sample = args.sample
        self.outdir = args.outdir
        self.thread = int(args.thread)
        self.bamdir = args.bam_dir
        self.well_sample = utils.two_col_to_dict(args.well_sample)
        self.index_dict = {}
        self.df = pd.DataFrame()

        # output files
        self.outstat = os.path.join(self.outdir, self.sample + ".substitution.csv")

    @utils.add_log
    def run(self):
        # check bam existence, skip wells with no data
        valid_samples = []
        bam_list = []
        for well, sample in self.well_sample.items():
            xbam = f"{self.bamdir}/{sample}/{sample}.PosTag.bam"
            if not os.path.exists(xbam):
                print(f"Well {well}: {sample} has no data, skipped.")
                continue
            valid_samples.append(sample)
            bam_list.append(xbam)

        if len(bam_list) == 0:
            raise ValueError("No valid bam files found for substitution.")

        mincpu = min(len(bam_list), self.thread)
        with Pool(mincpu) as pool:
            results = pool.map(Substitution.substitution_process, bam_list)

        for index, xout in zip(valid_samples, results):
            self.index_dict[index] = xout

        self.write_csv_file()
        self.sub_plot()
        self.add_help()

    @staticmethod
    def substitution_process(bam):
        for_base, rev_base, is_forward, is_reverse = Substitution.get_sub_tag(bam)
        xout = Substitution.sub_stat(for_base, rev_base, is_forward, is_reverse)
        return xout

    @staticmethod
    def get_sub_tag(bam):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, "rb")
        pysam.set_verbosity(save)
        is_reverse = {
            "cA": 0,
            "gA": 0,
            "tA": 0,
            "aC": 0,
            "gC": 0,
            "tC": 0,
            "aG": 0,
            "cG": 0,
            "tG": 0,
            "aT": 0,
            "cT": 0,
            "gT": 0,
        }
        is_forward = {
            "cA": 0,
            "gA": 0,
            "tA": 0,
            "aC": 0,
            "gC": 0,
            "tC": 0,
            "aG": 0,
            "cG": 0,
            "tG": 0,
            "aT": 0,
            "cT": 0,
            "gT": 0,
        }
        for_base = {"a": 0, "c": 0, "g": 0, "t": 0}
        rev_base = {"a": 0, "c": 0, "g": 0, "t": 0}
        snp_tags = [
            "",
            "cA",
            "gA",
            "tA",
            "aC",
            "gC",
            "tC",
            "aG",
            "cG",
            "tG",
            "aT",
            "cT",
            "gT",
        ]
        ref_tags = ["", "a", "c", "g", "t"]
        for read in bamfile.fetch(until_eof=True):
            try:
                snpmatch = re.match(
                    r"cA(\d+);gA(\d+);tA(\d+);aC(\d+);gC(\d+);tC(\d+);aG(\d+);cG(\d+);tG(\d+);aT(\d+);cT(\d+);gT(\d+);",
                    read.get_tag("SC"),
                    re.M,
                )
                totmatch = re.match(
                    r"a(\d+);c(\d+);g(\d+);t(\d+)", read.get_tag("TC"), re.M
                )
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

    @staticmethod
    def sub_stat(for_base, rev_base, is_forward, is_reverse):
        convertdict = {
            "a": ["aC", "aG", "aT"],
            "c": ["cA", "cG", "cT"],
            "g": ["gA", "gC", "gT"],
            "t": ["tA", "tC", "tG"],
        }
        subdict = {
            "a": "t",
            "t": "a",
            "c": "g",
            "g": "c",
            "aC": "tG",
            "aG": "tC",
            "aT": "tA",
            "cA": "gT",
            "cG": "gC",
            "cT": "gA",
            "gA": "cT",
            "gC": "cG",
            "gT": "cA",
            "tA": "aT",
            "tC": "aG",
            "tG": "aC",
        }
        outdict = {
            "aC": "A_to_C",
            "aG": "A_to_G",
            "aT": "A_to_T",
            "cA": "C_to_A",
            "cG": "C_to_G",
            "cT": "C_to_T",
            "gA": "G_to_A",
            "gC": "G_to_C",
            "gT": "G_to_T",
            "tA": "T_to_A",
            "tC": "T_to_C",
            "tG": "T_to_G",
        }
        outw = {}
        for x in ["a", "c", "g", "t"]:
            fbase = for_base[x]
            rbase = rev_base[subdict[x]]
            for y in convertdict[x]:
                fcov = is_forward[y] * 100 / float(fbase) if float(fbase) > 0 else 0
                rcov = (
                    is_reverse[subdict[y]] * 100 / float(rbase)
                    if float(rbase) > 0
                    else 0
                )
                outw[outdict[y]] = round(fcov + rcov, 3)

        return outw

    @utils.add_log
    def sub_plot(self):
        heatmap_df = self.df.T
        subplot = Substitution_plot_bulk(
            df_sub=heatmap_df, plot_type="heatmap"
        ).get_plotly_div()
        self.add_data(substitution_heatmap=subplot)

        row_sums = self.df.sum(axis=1)
        df_ratio = self.df.div(row_sums.replace(0, 1), axis=0)
        df_long = df_ratio.reset_index().melt(
            id_vars="index", var_name="type", value_name="ratio"
        )
        df_long = df_long.rename(columns={"index": "well"})
        subplot = Substitution_plot_bulk(
            df_sub=df_long, plot_type="bar", header=tc_header
        ).get_plotly_div()
        self.add_data(substitution_bar=subplot)

    @utils.add_log
    def write_csv_file(self):
        self.df = pd.DataFrame.from_dict(self.index_dict, orient="index")
        self.df = self.df[sorted_header]
        self.df.to_csv(self.outstat)

    @utils.add_log
    def add_help(self):
        self.add_help_content(
            name="substitution",
            content="Nucleotide substitution rate (plus strand + minus strand) per well and per type",
        )


@utils.add_log
def substitution(args):
    with Substitution(args) as runner:
        runner.run()


def get_opts_substitution(parser, sub_program):
    parser.add_argument(
        "--well_sample",
        help="tsv file of well numbers and sample names. The first column is well numbers and the second column is sample names.",
        required=True,
    )
    if sub_program:
        parser.add_argument(
            "--bam_dir", help="Required. bam dir from step conversion.", required=True
        )
        parser = s_common(parser)
    return parser
