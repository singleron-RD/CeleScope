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


class SMK:

    def __init__(self, SMK_read2, match_barcode_file, sample, SMK_barcode_fasta, tsne_file, UMI_min, percent_min, 
        combine_cluster, run_summary, UMI2, outdir):
        self.SMK_read2 = SMK_read2
        self.match_barcode_file = match_barcode_file
        self.sample = sample
        self.SMK_barcode_fasta = SMK_barcode_fasta
        self.tsne_file = tsne_file
        self.percent_min = percent_min
        self.UMI_min = UMI_min
        self.combine_cluster = combine_cluster
        self.run_summary = run_summary
        self.UMI2 = UMI2
        self.outdir = outdir

        self.SMK_barcode_dic = {}
        self.SMK_barcode_length = None
        self.res_dic = SMK.genDict()
        self.res_sum_dic = SMK.genDict(dim=2)
        self.cell_total = 0
        self.cell_assigned = 0
        self.cell_failed = 0
        self.cell_ambiguity = 0
        self.reads = 0
        self.reads_unmapped = 0
        self.reads_mapped_not_in_cell = 0
        self.reads_mapped_in_cell = 0
        self.reads_too_short = 0
        self.match_barcode = []
        self.out_df = f'{outdir}/{sample}_UMI.tsv'
        self.res_df_file = f'{outdir}/{sample}_read_count.tsv'
        self.out_UMI2_df = f'{outdir}/{sample}_UMI2.tsv'
        self.df_tsne_summary_file = f'{outdir}/{sample}_tag_summary.tsv'
        self.df_count_file = f'{outdir}/{sample}_cluster_count.tsv'
        self.cluster_plot = f'{outdir}/{sample}_cluster_plot.pdf'
        if combine_cluster:
            self.df_count_combine_file = f'{outdir}/{sample}_combine_cluster_count.tsv'
            self.combine_cluster_plot = f'{outdir}/{sample}_combine_cluster_plot.pdf'

        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)

    @staticmethod
    def genDict(dim=3):
        if dim==1:
            return defaultdict(int)
        else:
            return defaultdict(lambda: SMK.genDict(dim-1))

    def read_barcode_file(self):
        with open(self.match_barcode_file, "rt") as f:
            while True:
                line = f.readline().strip()
                if not line:
                    break
                self.match_barcode.append(line)    
        self.cell_total = len(self.match_barcode)          

    @staticmethod
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
                if SMK.hamming_distance(seq, SMK_barcode_seq) < 3:
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
        res_df = pd.DataFrame(rows)
        res_df.rename(columns={0:"barcode", 1:"SMK_barcode_name", 2:"UMI", 3:"read_count"}, inplace=True)
        res_df.to_csv(self.res_df_file, sep="\t", index=False)

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
        stat_file = self.outdir + "/stat.txt"
        stat.to_csv(stat_file, sep=":", header=None, index=False)
        
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

    @staticmethod
    def tag_type(sub_df, UMI_min, percent_min):
        UMI_sum = sub_df["UMI_count"].sum()
        UMI_max = sub_df["UMI_count"].max()
        UMI_max_percent = UMI_max/UMI_sum
        if UMI_sum < UMI_min:
            return "Undetermined"
        elif  UMI_max_percent < percent_min:
            return "Multiplet"
        else:
            return sub_df.loc[sub_df["UMI_count"] == UMI_max, "SMK_barcode"].values[0]

    @staticmethod
    def write_and_plot(df, column_name, count_file, plot_file):
        df_count = df.groupby(["tag",column_name]).size().unstack()
        df_count.fillna(0, inplace=True)
        df_count.to_csv(count_file, sep="\t")
        df_percent = df_count/df_count.sum()
        df_plot = df_percent.stack().reset_index()
        df_plot.rename({0:"percent"}, axis=1, inplace=True)

        #plot
        colors = list(mpl.colors.cnames.keys())
        fig, ax = plt.subplots(figsize=(20, 10))  
        types = df_plot["tag"].drop_duplicates()
        margin_bottom = np.zeros(len(df_plot[column_name].drop_duplicates()))

        for num, type in enumerate(types):
            values = list(df_plot.loc[df_plot["tag"]==type, "percent"])
            df_plot[df_plot['tag'] == type].plot.bar(x=column_name, y='percent', ax=ax, stacked=True, 
                                            bottom = margin_bottom, label=type,color=colors[num*3+1])
            margin_bottom += values
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.title("SMK tag fraction")
        fig.savefig(plot_file)

    def tag_count(self):
        if not self.UMI2:
            for barcode in self.res_dic:
                for SMK_barcode_name in self.res_dic[barcode]:
                    self.res_sum_dic[barcode][SMK_barcode_name] = len(self.res_dic[barcode][SMK_barcode_name])
        else:
            for barcode in self.res_dic:
                for SMK_barcode_name in self.res_dic[barcode]:
                    UMI2_count = sum(list(map(lambda x:x > 1, list(self.res_dic[barcode][SMK_barcode_name].values()))))
                    self.res_sum_dic[barcode][SMK_barcode_name] = UMI2_count                    
        self.df = pd.DataFrame(self.res_sum_dic)
        self.df = self.df.T
        self.df.fillna(0, inplace=True)
        self.df.to_csv(self.out_df, sep="\t")

    def summary(self):

        if not self.run_summary:
            self.tag_count()
        else:
            if not self.UMI2:
                df_file = glob.glob(sample + "_UMI.tsv")[0]
                self.df = pd.read_csv(df_file, sep="\t", index_col=0)
                self.df.fillna(0, inplace=True)
            else:
                df_read_file = glob.glob(sample+"_read_count.tsv")[0]
                df_read = pd.read_csv(df_read_file, sep="\t", index_col=0)
                df_read_count = df_read[["SMK_barcode_name", "read_count"]]
                df_UMI2_count = df_read_count[df_read_count["read_count"] > 1].reset_index()
                df_UMI2_count = df_UMI2_count.groupby(["barcode", "SMK_barcode_name"]).agg("count")
                df_UMI2_count.rename(columns={"read_count": "UMI_count"}, inplace=True)
                self.df = df_UMI2_count.unstack("SMK_barcode_name", fill_value=0)
                self.df.columns = self.df.columns.droplevel()
                self.df.fillna(0,inplace=True)
                self.df.to_csv(self.out_UMI2_df, sep="\t")

        df_tsne = pd.read_csv(self.tsne_file, sep="\t", index_col=0)

        df_tag_UMI = self.df.stack().reset_index()
        df_tag_UMI.columns = ["barcode", "SMK_barcode", "UMI_count"]

        barcode_tag_type = df_tag_UMI.groupby(["barcode"]).apply(SMK.tag_type, UMI_min=self.UMI_min,
            percent_min=self.percent_min)
        df_summary = pd.DataFrame(barcode_tag_type)
        df_summary.columns = ["tag"]
        df_tsne_summary = pd.merge(df_tsne, df_summary, how="left", left_index=True, right_index=True)
        df_tsne_summary.fillna({"tag": "Undetermined"}, inplace=True)
        if self.combine_cluster:
            df_combine_cluster = pd.read_csv(self.combine_cluster,sep="\t", header=None)
            df_combine_cluster.columns = ["cluster", "combine_cluster"]
            df_tsne_combine_cluster_summary = pd.merge(df_tsne_summary, df_combine_cluster,
                on=["cluster"], how="left", left_index=True).set_index(df_tsne_summary.index)
            df_tsne_combine_cluster_summary.to_csv(self.df_tsne_summary_file, sep="\t")
        else:
            df_tsne_summary.to_csv(self.df_tsne_summary_file, sep="\t")
        
        SMK.write_and_plot(df=df_tsne_summary, column_name="cluster", count_file=self.df_count_file,
            plot_file=self.cluster_plot)
        
        if self.combine_cluster:
            SMK.write_and_plot(df=df_tsne_combine_cluster_summary, column_name="combine_cluster",
            count_file=self.df_count_combine_file, plot_file=self.combine_cluster_plot)
        
    
    def run(self):
        logger1.info("running")
        if not self.run_summary:
            self.read_SMK_barcode()
            self.read_barcode_file()
            self.read_to_dic()
        self.summary()
        logger1.info("done")


def get_opts_demultiplex(parser, sub_program):
    parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path")
    parser.add_argument("--SMK_barcode_fasta", help="SMK barcode fasta")
    parser.add_argument("--UMI_min", help="keep tags have UMI>=UMI_min", default = 50)
    parser.add_argument("--percent_min", help="keep tags have percent>=percent_min", default = 0.7)
    parser.add_argument("--combine_cluster", help="conbine cluster tsv file", default = None)
    parser.add_argument("--UMI2", help="only consider UMI with as least two read support", action='store_true')
    parser.add_argument("--run_summary", help="only run demultiplex summary", action='store_true')
    if sub_program:
        parser.add_argument("--SMK_read2", help="SMK clean read2", required=True)
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)


def demultiplex(args):
    match_dir = args.match_dir
    sample = args.sample
    outdir = args.outdir
    SMK_barcode_fasta = args.SMK_barcode_fasta
    UMI_min = int(args.UMI_min)
    percent_min = float(args.percent_min)
    combine_cluster = args.combine_cluster
    run_summary = args.run_summary
    UMI2 = args.UMI2
    SMK_read2 = args.SMK_read2

    match_barcode_file1 = glob.glob("{match_dir}/05.count/*_cellbarcode.tsv".format(match_dir=match_dir))
    match_barcode_file2 = glob.glob("{match_dir}/05.count/matrix_10X/*_cellbarcode.tsv".format(match_dir=match_dir))
    match_barcode_file = (match_barcode_file1 + match_barcode_file2)[0]
    tsne_file = glob.glob("{match_dir}/06.analysis/tsne_coord.tsv".format(match_dir=match_dir))[0]

    s = SMK(SMK_read2=SMK_read2, match_barcode_file=match_barcode_file, sample=sample, 
        SMK_barcode_fasta=SMK_barcode_fasta, tsne_file=tsne_file, UMI_min=UMI_min, percent_min=percent_min, 
        combine_cluster=combine_cluster, run_summary=run_summary, UMI2=UMI2, outdir=outdir)
    s.run()