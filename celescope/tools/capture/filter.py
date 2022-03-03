import json

import pandas as pd

from celescope.tools import utils
from celescope.tools.count import Count
from celescope.tools.step import Step, s_common
from celescope.tools.capture.threshold import Threshold
from celescope.__init__ import HELP_DICT
from celescope.tools.capture.__init__ import SUM_UMI_COLNAME
from celescope.tools.analysis_wrapper import Report_runner


def get_opts_filter(parser, sub_program):
    parser.add_argument('--not_correct_UMI', help='Do not perform UMI correction.', action='store_true')

    parser.add_argument(
        "--read_threshold_method",
        help='method to find read threshold. UMIs with `support reads` < `read threshold` are filtered.',
        choices=['otsu', 'auto', 'hard', 'none'],
        default='otsu'
    )
    parser.add_argument(
        "--read_hard_threshold",
        help='int, use together with `--read_threshold_method hard`',
    )

    parser.add_argument(
        "--umi_threshold_method",
        help='method to find UMI threshold. Cell barcode with `UMI` < `UMI threshold` are considered negative.',
        choices=['otsu', 'auto', 'hard', 'none'],
        default='auto'
    )
    parser.add_argument(
        "--umi_hard_threshold",
        help='int, use together with `--umi_threshold_method hard`',
    )

    if sub_program:
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--raw_read_count_file', required=True)
        s_common(parser)


class Filter(Step):
    def __init__(self, args, display_title='Filtering'):
        super().__init__(args, display_title)

        # data
        with open(args.raw_read_count_file) as f:
            self.count_dict = json.load(f)

        self.raw_umi = 0
        self.total_corrected_umi = 0
        self.del_umi = 0
        self.read_threshold_dict = {}
        self.umi_threshold_dict = {}  # if not set explicitly, use 1 as default

        self.barcode_ref_umi_dict = utils.genDict(dim=2)
        self.ref_barcode_umi_dict = utils.genDict(dim=2)

        report_runner = Report_runner(args)
        df_tsne, _df_marker = report_runner.get_df()
        self.df_filter_tsne = df_tsne.copy()

        # out
        self.corrected_read_count_file = f'{self.out_prefix}_corrected_read_count.json'
        self.filter_read_count_file = f'{self.out_prefix}_filtered_read_count.json'
        self.filter_tsne_file = f'{self.out_prefix}_filtered_UMI_tsne.csv'

    @utils.add_log
    def correct_umi(self):
        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                self.raw_umi += len(self.count_dict[barcode][ref])
                n_corrected_umi, _n_corrected_read = Count.correct_umi(self.count_dict[barcode][ref])
                if self.debug:
                    print(f'{barcode} {ref} {n_corrected_umi}')
                self.total_corrected_umi += n_corrected_umi

        self.add_metric(
            name='Number of Raw UMI',
            value=self.raw_umi,
            help_info='number of total raw UMI',
        )

        self.add_metric(
            name='Number of Corrected UMI',
            value=self.total_corrected_umi,
            total=self.raw_umi,
            help_info='correct sequencing errors in the UMI sequences ',
        )

    @utils.add_log
    def write_correct_umi_json(self):
        with open(self.corrected_read_count_file, 'w') as fp:
            json.dump(self.count_dict, fp, indent=4)

    @utils.add_log
    def get_read_threshold(self):
        self.add_metric('Read Threshold Method', self.args.read_threshold_method)

        read_dict = utils.genDict(dim=1, valType=list)

        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                read_dict[ref] += list(self.count_dict[barcode][ref].values())

        if self.debug:
            print(read_dict)

        for ref in read_dict:              
            otsu_plot_path = f'{self.out_prefix}_{ref}_read_otsu.png'
            runner = Threshold(
                array=read_dict[ref], 
                threshold_method=self.args.read_threshold_method, 
                otsu_plot_path=otsu_plot_path,
                hard_threshold=self.args.read_hard_threshold
            )
            read_threshold = runner.run()
            
            self.read_threshold_dict[ref] = read_threshold
            self.add_metric(
                f'{ref} Read Threshold', 
                read_threshold,
                help_info='filter UMI with less than this number of reads',
            )


    @utils.add_log
    def filter_read(self):
        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                for umi in self.count_dict[barcode][ref]:
                    read_count = self.count_dict[barcode][ref][umi]
                    if read_count < self.read_threshold_dict[ref]:
                        self.del_umi += 1
                        self.count_dict[barcode][ref][umi] = 0

        """
        self.add_metric(
            name='Number of Total Filtered UMI',
            value=self.del_umi,
            total=self.raw_umi,
        )
        """

    def write_filter_read_json(self):
        with open(self.filter_read_count_file, 'w') as fp:
            json.dump(self.count_dict, fp, indent=4)


    @utils.add_log
    def set_barcode_ref_umi_dict(self):

        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                for umi in self.count_dict[barcode][ref]:
                    if self.count_dict[barcode][ref][umi] > 0:
                        self.barcode_ref_umi_dict[barcode][ref] += 1
        if self.debug:
            print(dict(self.barcode_ref_umi_dict))

    @utils.add_log
    def set_ref_barcode_umi_dict(self):

        for barcode in self.barcode_ref_umi_dict:
            for ref in self.barcode_ref_umi_dict[barcode]:
                self.ref_barcode_umi_dict[ref][barcode] = self.barcode_ref_umi_dict[barcode][ref]

    def get_umi_threshold(self):
        self.add_metric('UMI Threshold Method', self.args.umi_threshold_method)

        for ref in self.ref_barcode_umi_dict:
            umi_array = list(self.ref_barcode_umi_dict[ref].values())
            otsu_plot_path = f'{self.out_prefix}_{ref}_UMI_otsu.png'
            runner = Threshold(
                umi_array, 
                threshold_method=self.args.umi_threshold_method, 
                otsu_plot_path=otsu_plot_path,
                hard_threshold=self.args.umi_hard_threshold
            )
            umi_threshold = runner.run()
            
            umi_threshold = max(1, umi_threshold)
            self.umi_threshold_dict[ref] = umi_threshold
            self.add_metric(f'{ref} UMI Threshold', umi_threshold)


    @utils.add_log
    def filter_umi(self):
        for ref in self.ref_barcode_umi_dict:
            for barcode in self.ref_barcode_umi_dict[ref]:
                if self.ref_barcode_umi_dict[ref][barcode] < self.umi_threshold_dict[ref]:
                    self.ref_barcode_umi_dict[ref][barcode] = 0

    @utils.add_log
    def add_umi_write_tsne(self):
        for ref in self.umi_threshold_dict:
            self.df_filter_tsne[ref] = pd.Series(self.ref_barcode_umi_dict[ref])
            self.df_filter_tsne[ref].fillna(0, inplace=True)

        refs = list(self.umi_threshold_dict.keys()) 
        self.df_filter_tsne[SUM_UMI_COLNAME] = self.df_filter_tsne[refs].sum(axis=1)
        self.df_filter_tsne.to_csv(self.filter_tsne_file)


    def add_some_metrics(self):
        df = self.df_filter_tsne

        cell_total = len(df)
        df_positive = df[df[SUM_UMI_COLNAME] > 0]
        n_cell_positive = len(df_positive)
        self.add_metric(
            name='Number of Positive Cells after Filtering',
            value=n_cell_positive,
            total=cell_total,
        )

    def run(self):
        if not self.args.not_correct_UMI:
            self.correct_umi()
            self.write_correct_umi_json()

        self.get_read_threshold()
        self.filter_read()
        self.write_filter_read_json()

        self.set_barcode_ref_umi_dict()
        self.set_ref_barcode_umi_dict()
        self.get_umi_threshold()
        self.filter_umi()

        self.add_umi_write_tsne()
        self.add_some_metrics()
