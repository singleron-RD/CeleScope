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
from celescope.tools.utils import format_number
from celescope.tools.report import reporter

logger1 = logging.getLogger(__name__)

def genDict(dim=3):
    if dim==1:
        return defaultdict(int)
    else:
        return defaultdict(lambda: genDict(dim-1))


def hamming_distance(string1, string2):
    distance = 0
    length = len(string1)
    length2 = len(string2)
    if (length != length2):
        raise Exception("string1 and string2 do not have same length")
    for i in range(length):
        if string1[i] != string2[i]:
            distance += 1
    return distance


class smk_mapping:

    def __init__(self, SMK_read2, match_barcode_file, sample, SMK_barcode_fasta, outdir):
        self.SMK_read2 = SMK_read2
        self.match_barcode_file = match_barcode_file
        self.sample = sample
        self.SMK_barcode_fasta = SMK_barcode_fasta
        self.outdir = outdir

        self.SMK_barcode_dic = {}
        self.SMK_barcode_length = None

        self.res_dic = genDict()
        self.res_sum_dic = genDict(dim=2)

        self.reads = 0
        self.reads_unmapped = 0
        self.reads_mapped_not_in_cell = 0
        self.reads_mapped_in_cell = 0
        self.reads_too_short = 0
        self.match_barcode = []
        self.cell_read_count_file = f'{outdir}/{sample}_cell_read_count.tsv'
        self.cell_umi_count_file = f'{outdir}/{sample}_cell_UMI_count.tsv'
        self.stat_file = self.outdir + "/stat.txt"

        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)

    def read_barcode_file(self):
        with open(self.match_barcode_file, "rt") as f:
            while True:
                line = f.readline().strip()
                if not line:
                    break
                self.match_barcode.append(line)    
        self.cell_total = len(self.match_barcode)

    def read_to_dic(self):

        read2 = gzip.open(self.SMK_read2, "rt")
        index = 0
        seq_length = len(list(self.SMK_barcode_dic.values())[0])
        while True:
            line1 = read2.readline()
            line2 = read2.readline()
            line3 = read2.readline()
            line4 = read2.readline()
            if not line4:
                break
            self.reads += 1
            if self.reads % 1000000 == 0:
                logger1.info(str(self.reads) + " reads done.")
            attr = str(line1).strip("@").split("_")
            barcode = str(attr[0])
            umi = str(attr[1])
            seq = line2.strip()
            if len(seq) < self.SMK_barcode_length:
                self.reads_too_short += 1
                self.reads_unmapped += 1
                continue
            seq = seq[0:self.SMK_barcode_length]      
            mapped = False
            for SMK_barcode_name in self.SMK_barcode_dic:
                SMK_barcode_seq = self.SMK_barcode_dic[SMK_barcode_name]
                if hamming_distance(seq, SMK_barcode_seq) < 3:
                    mapped = True
                    if not (barcode in self.match_barcode):
                        self.reads_mapped_not_in_cell += 1
                        continue                    
                    self.res_dic[barcode][SMK_barcode_name][umi] += 1
                    self.reads_mapped_in_cell += 1
                    continue
            if not mapped:
                self.reads_unmapped += 1
        
        # write dic to pandas df
        rows = []
        for barcode in self.res_dic:
            for SMK_barcode_name in self.res_dic[barcode]:
                for umi in self.res_dic[barcode][SMK_barcode_name]:
                    rows.append([barcode, SMK_barcode_name, umi, self.res_dic[barcode][SMK_barcode_name][umi]])
        df_cell_read_count = pd.DataFrame(rows)
        df_cell_read_count.rename(columns={0:"barcode", 1:"SMK_barcode_name", 2:"UMI", 3:"read_count"}, inplace=True)
        df_cell_read_count.to_csv(self.cell_read_count_file, sep="\t", index=False)

        # stat
        stat = pd.DataFrame(
        {"item": [
            "Reads Mapped in Cell", 
            "Reads Mapped not in Cell", 
            "Reads Unmapped",
            "Reads Unmapped too Short", 
            ],
        "count": [
            self.reads_mapped_in_cell,
            self.reads_mapped_not_in_cell,
            self.reads_unmapped,
            self.reads_too_short,
            ],
        }, columns=["item", "count"])
        stat["count_percent"] = stat["count"].apply(lambda x:f'{x}({round(x/self.reads * 100, 2)}%)')
        stat = stat.loc[:, ["item", "count_percent"]]
        stat.to_csv(self.stat_file, sep=":", header=None, index=False)
        
    def read_SMK_barcode(self):

        with open(self.SMK_barcode_fasta, "rt") as f:
            while True:
                line1 = f.readline()
                line2 = f.readline()
                if not line2:
                    break
                barcode_name = line1.strip(">").strip()
                barcode_seq = line2.strip()
                SMK_barcode_length = len(barcode_seq)
                if not self.SMK_barcode_length:
                    self.SMK_barcode_length = SMK_barcode_length
                else:
                    if self.SMK_barcode_length != SMK_barcode_length:
                        raise Exception("SMK barcode have different length")
                self.SMK_barcode_dic[barcode_name] = barcode_seq

    def tag_count(self):
        for barcode in self.res_dic:
                for SMK_barcode_name in self.res_dic[barcode]:
                    self.res_sum_dic[barcode][SMK_barcode_name] = len(self.res_dic[barcode][SMK_barcode_name])
                
        df_cell_umi_count = pd.DataFrame(self.res_sum_dic)
        df_cell_umi_count = df_cell_umi_count.T
        df_cell_umi_count.fillna(0, inplace=True)
        df_cell_umi_count = df_cell_umi_count.astype(int)
        df_cell_umi_count.to_csv(self.cell_umi_count_file, sep="\t")

    def run(self):
        logger1.info('mapping smk start...!')
        self.read_SMK_barcode()
        self.read_barcode_file()
        self.read_to_dic()
        self.tag_count()
        t = reporter(name='mapping_smk', assay="smk", sample=self.sample, 
        stat_file=self.stat_file, outdir=self.outdir + '/..')
        t.get_report()
        logger1.info('mapping smk done...!')


def get_opts_mapping_smk(parser, sub_program):
    parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path")
    parser.add_argument("--SMK_barcode_fasta", help="SMK barcode fasta")
    if sub_program:
        parser.add_argument("--SMK_read2", help="SMK clean read2", required=True)
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)


def mapping_smk(args):
    match_dir = args.match_dir
    sample = args.sample
    outdir = args.outdir
    SMK_barcode_fasta = args.SMK_barcode_fasta
    SMK_read2 = args.SMK_read2

    match_barcode_file1 = glob.glob("{match_dir}/05.count/*_cellbarcode.tsv".format(match_dir=match_dir))
    match_barcode_file2 = glob.glob("{match_dir}/05.count/matrix_10X/*_cellbarcode.tsv".format(match_dir=match_dir))
    match_barcode_file = (match_barcode_file1 + match_barcode_file2)[0]

    s = smk_mapping(SMK_read2, match_barcode_file, sample, SMK_barcode_fasta, outdir)
    s.run()
