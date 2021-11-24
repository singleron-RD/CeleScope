import glob
import os

import numpy as np
import pandas as pd

from celescope.tools.report import reporter
from celescope.tools.utils import add_log, format_stat, read_barcode_file


class Count_cite():

    def __init__(
        self,
        sample,
        outdir,
        assay,
        read_count_file,
        match_dir,
    ):
        self.sample = sample
        self.outdir = outdir
        self.assay = assay
        self.read_count_file = read_count_file
        self.match_dir = match_dir
        self.match_barcode, self.cell_total = read_barcode_file(match_dir)
        self.df_read_count = pd.read_csv(read_count_file, sep="\t", index_col=0)
        self.tsne_file = glob.glob(f'{match_dir}/*analysis/*tsne_coord.tsv')[0]

        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)

        # out
        self.mtx = f'{outdir}/{sample}_citeseq.mtx.gz'
        self.stats = None
        self.stat_file = f'{self.outdir}/stat.txt'

    @add_log
    def run(self):
        stats = pd.Series()
        mapped_read = self.df_read_count['read_count'].sum()

        # in cell
        df_read_count_in_cell = self.df_read_count[self.df_read_count.index.isin(self.match_barcode)]
        mapped_read_in_cell = int(df_read_count_in_cell['read_count'].sum())
        stats = stats.append(pd.Series(
            format_stat(mapped_read_in_cell, mapped_read),
            index=['Mapped Reads in Cells']
        ))

        # UMI
        df_UMI_in_cell = df_read_count_in_cell.reset_index().groupby([
            'barcode', 'barcode_name']).agg({'UMI': 'count'})
        df_UMI_in_cell = df_UMI_in_cell.reset_index()
        df_UMI_in_cell = df_UMI_in_cell.pivot(
            index='barcode', columns='barcode_name', values='UMI')
        df_cell = pd.DataFrame(index=self.match_barcode)
        df_UMI_cell = pd.merge(
            df_cell,
            df_UMI_in_cell,
            how="left",
            left_index=True,
            right_index=True)

        # fillna
        df_UMI_cell.fillna(0, inplace=True)
        df_UMI_cell = df_UMI_cell.astype(int)
        df_UMI_cell_out = df_UMI_cell.T
        df_UMI_cell_out.to_csv(self.mtx, sep='\t', compression='gzip')

        # UMI
        UMIs = df_UMI_cell.apply(sum, axis=1)
        median = round(np.median(UMIs), 2)
        mean = round(np.mean(UMIs), 2)
        stats = stats.append(pd.Series(
            str(median),
            index=['Median UMI per Cell']
        ))

        stats = stats.append(pd.Series(
            str(mean),
            index=['Mean UMI per Cell']
        ))

        self.stats = stats

    def report(self):

        self.stats.to_csv(self.stat_file, sep=':', header=False)
        t = reporter(
            name='count_cite',
            assay=self.assay,
            sample=self.sample,
            stat_file=self.stat_file,
            outdir=self.outdir + '/..')
        t.get_report()

def get_opts_count_cite(parser, sub_program):
    parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path", required=True)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument("--read_count_file", help="tag read count file")


def count_cite(args):

    count_cite_object = Count_cite(
        args.sample,
        args.outdir,
        args.assay,
        args.read_count_file,
        args.match_dir,
    )
    count_cite_object.run()
    count_cite_object.report()
