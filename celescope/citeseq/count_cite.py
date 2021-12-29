import numpy as np
import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT


class Count_cite(Step):

    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        self.df_read_count = pd.read_csv(args.read_count_file, sep='\t', index_col=0)

        self.tsne_dict = utils.parse_match_dir(args.match_dir)
        self.match_barcode = self.tsne_dict['match_barcode']

        # out
        self.mtx = f'{self.out_prefix}_citeseq.mtx.gz'

    @utils.add_log
    def run(self):
        mapped_read = int(self.df_read_count['read_count'].sum())

        # in cell
        df_read_count_in_cell = self.df_read_count[self.df_read_count.index.isin(self.match_barcode)]
        mapped_read_in_cell = int(df_read_count_in_cell['read_count'].sum())
        self.add_metric(
            name='Mapped Reads in Cells',
            value=mapped_read_in_cell,
            total=mapped_read,
        )

        # UMI
        df_UMI_in_cell = df_read_count_in_cell.reset_index().groupby([
            'barcode', 'tag_name']).agg({'UMI': 'count'})
        df_UMI_in_cell = df_UMI_in_cell.reset_index()
        df_UMI_in_cell = df_UMI_in_cell.pivot(
            index='barcode', columns='tag_name', values='UMI')
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
        median_umi = round(np.median(UMIs), 2)
        mean_umi = round(np.mean(UMIs), 2)
        self.add_metric(
            name='Median UMI per Cell',
            value=float(median_umi),
        )
        self.add_metric(
            name='Mean UMI per Cell',
            value=float(mean_umi),
        )


def get_opts_count_cite(parser, sub_program):
    if sub_program:
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'], required=True)
        parser.add_argument("--read_count_file", help="tag read count file")
        s_common(parser)


def count_cite(args):

    with Count_cite(args, display_title='Count') as runner:
        runner.run()
