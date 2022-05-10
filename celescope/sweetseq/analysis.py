import json

import pandas as pd
import numpy as np

from celescope.tools import analysis_wrapper
from celescope.tools.plotly_plot import Tsne_plot
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT


class Analysis(Step):
    """
    ## Features


    ## Output

    """

    def __init__(self, args, display_title=None):

        super().__init__(args, display_title)
     
        # input
        with open(args.raw_read_count_file) as f:
            self.read_count_dict = json.load(f)

        self.display_title = display_title
        self.ref_barcode_umi_dict = utils.genDict(dim=2)

        report_runner = analysis_wrapper.Report_runner(self.args, display_title=self.display_title)
        self.df_tsne, self.df_marker = report_runner.get_df()
        self.df_tsne_umi = self.df_tsne.copy()

        # data
        self.ref_name = None

        # const 
        self._table_id = 'marker_genes'

        # out
        self.tsne_umi_file = f'{self.out_prefix}_UMI_tsne.csv'


    @utils.add_log
    def set_ref_barcode_umi_dict(self):

        for barcode in self.read_count_dict:
            for ref in self.read_count_dict[barcode]:
                for umi in self.read_count_dict[barcode][ref]:
                    if self.read_count_dict[barcode][ref][umi] > 0:
                        self.ref_barcode_umi_dict[ref][barcode] += 1

    @utils.add_log
    def add_umi_write_tsne(self):
        for ref in self.ref_barcode_umi_dict:
            self.df_tsne_umi[ref] = pd.Series(self.ref_barcode_umi_dict[ref])
            self.df_tsne_umi[ref].fillna(0, inplace=True)
            self.df_tsne_umi[ref + '_log2'] = np.log2(self.df_tsne_umi[ref] + 1)

        self.df_tsne_umi.to_csv(self.tsne_umi_file)

        if len(self.ref_barcode_umi_dict) > 1:
            raise ValueError("Error: More than one barcode in tag_barcode_fasta. Currently, Multiple barcodes are not supported to display in HTML report.")
        else:
            self.ref_name = list(self.ref_barcode_umi_dict.keys())[0] + '_log2'

   
    def add_html_data(self):
        tsne_cluster = Tsne_plot(self.df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_plot = Tsne_plot(self.df_tsne_umi, self.ref_name, discrete=False)
        # tsne_plot.set_color_scale(['LightGrey', 'Red'])
        tsne_umi = tsne_plot.get_plotly_div()
        self.add_data(tsne_umi=tsne_umi)

        table_dict = self.get_table_dict(
            title='Marker Genes by Cluster',
            table_id=self._table_id,
            df_table=self.df_marker,
        )
        self.add_data(table_dict=table_dict)

    def run(self):
        self.set_ref_barcode_umi_dict()
        self.add_umi_write_tsne()
        self.add_html_data()


@utils.add_log
def analysis(args):
    with Analysis(args, display_title='Analysis') as runner:
        runner.run()


def get_opts_analysis(parser, sub_program):
    if sub_program:
        parser.add_argument('--raw_read_count_file', help='raw read', required=True)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        s_common(parser)