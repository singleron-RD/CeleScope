from celescope.tools.report import reporter
from celescope.tools.utils import format_number, log
from celescope.tools.utils import read_fasta, process_read, gen_stat
from celescope.tools.barcode import parse_pattern
import gzip
import os
import pandas as pd
import logging
import numpy as np
import sys
import argparse
import re
import glob
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def genDict(dim=3):
    if dim == 1:
        return defaultdict(int)
    else:
        return defaultdict(lambda: genDict(dim - 1))


class smk_mapping:

    def __init__(
        self,
        SMK_read2,
        sample,
        SMK_pattern,
        SMK_linker,
        SMK_barcode,
        outdir,):

        # read args
        self.SMK_read2 = SMK_read2
        self.sample = sample
        self.SMK_pattern = SMK_pattern
        self.SMK_barcode = SMK_barcode
        self.SMK_linker = SMK_linker
        self.outdir = outdir

        # process
        self.SMK_barcode_dic, self.SMK_barcode_length = read_fasta(self.SMK_barcode, equal=True)
        self.SMK_linker_dic, self.SMK_linker_length = read_fasta(self.SMK_linker, equal=True)
        self.pattern_dict = parse_pattern(self.SMK_pattern)

        # check barcode length
        barcode1 = self.pattern_dict["C"][0]
        pattern_barcode_length = barcode1[1] - barcode1[0]
        if pattern_barcode_length != self.SMK_barcode_length:
            raise Exception(
                f'''SMK barcode fasta length {self.SMK_barcode_length} 
                != pattern barcode length {pattern_barcode_length}'''
            )

        self.res_dic = genDict()
        self.res_sum_dic = genDict(dim=2)

        self.match_barcode = []
        self.read_count_file = f'{outdir}/{sample}_read_count.tsv'
        self.UMI_count_file = f'{outdir}/{sample}_UMI_count.tsv'
        self.stat_file = self.outdir + "/stat.txt"

        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)

    def read_to_dic(self):

        self.res_dic, self.metrics = process_read(
            self.SMK_read2,
            self.pattern_dict,
            self.SMK_barcode_dic,
            self.SMK_linker_dic,
            self.SMK_barcode_length,
            self.SMK_linker_length,
        )

        # write dic to pandas df
        rows = []
        for barcode in self.res_dic:
            for SMK_barcode_name in self.res_dic[barcode]:
                for umi in self.res_dic[barcode][SMK_barcode_name]:
                    rows.append([barcode, SMK_barcode_name, umi,
                                 self.res_dic[barcode][SMK_barcode_name][umi]])
        df_read_count = pd.DataFrame(rows)
        df_read_count.rename(
            columns={
                0: "barcode",
                1: "SMK_barcode_name",
                2: "UMI",
                3: "read_count"},
            inplace=True)
        df_read_count.to_csv(
            self.read_count_file, sep="\t", index=False)

        # stat
        stat = pd.DataFrame(
            {
                "item": [
                    "Reads Mapped",
                    "Reads Unmapped too Short",
                    "Reads Unmapped Invalid Linker",
                    "Reads Unmapped Invalid Barcode",
                ],
                "count": [
                    self.metrics["Reads Mapped"],
                    self.metrics["Reads Unmapped too Short"],
                    self.metrics["Reads Unmapped Invalid Linker"],
                    self.metrics["Reads Unmapped Invalid Barcode"],
                ],
            }, 
            columns=["item", "count"]
        )
        stat['total_count'] = self.metrics["Total Reads"]
        gen_stat(stat, self.stat_file)

    def tag_count(self):
        for barcode in self.res_dic:
            for SMK_barcode_name in self.res_dic[barcode]:
                self.res_sum_dic[barcode][SMK_barcode_name] = len(
                    self.res_dic[barcode][SMK_barcode_name])

        df_umi_count = pd.DataFrame(self.res_sum_dic)
        df_umi_count = df_umi_count.T
        df_umi_count.fillna(0, inplace=True)
        df_umi_count = df_umi_count.astype(int)
        df_umi_count.to_csv(self.UMI_count_file, sep="\t")

    @log
    def run(self):
        self.read_to_dic()
        self.tag_count()
        t = reporter(
            name='mapping_smk', assay="smk", sample=self.sample,
            stat_file=self.stat_file, outdir=self.outdir + '/..')
        t.get_report()


def get_opts_mapping_smk(parser, sub_program):
    parser.add_argument("--SMK_pattern", help="SMK read2 pattern")
    parser.add_argument("--SMK_linker", help="SMK read2 linker fasta path")
    parser.add_argument("--SMK_barcode", help="SMK read2 barcode fasta path ")
    if sub_program:
        parser.add_argument(
            "--SMK_read2",
            help="SMK clean read2",
            required=True)
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)


def mapping_smk(args):
    sample = args.sample
    outdir = args.outdir
    SMK_pattern = args.SMK_pattern
    SMK_linker = args.SMK_linker
    SMK_barcode = args.SMK_barcode
    SMK_read2 = args.SMK_read2

    s = smk_mapping(
        SMK_read2,
        sample,
        SMK_pattern,
        SMK_linker,
        SMK_barcode,
        outdir)
    s.run()
