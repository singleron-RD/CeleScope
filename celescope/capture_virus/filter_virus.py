import json
from collections import defaultdict

import numpy as np
import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.count import Count
from celescope.tools.step import Step, s_common
import celescope.capture_virus.otsu as otsu 
from celescope.__init__ import HELP_DICT


def get_opts_filter_virus(parser, sub_program):

    parser.add_argument('--min_support_reads', help='Minimum number of reads to support a UMI', default=2)

    parser.add_argument(
        "--umi_threshold_method", 
        help='method to find virus UMI threshold',
        choices=['otsu', 'auto', 'hard'], 
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

def filter_virus(args):
    
    with Filter_virus(args, display_title='Filtering') as runner:
        runner.run()


class Filter_virus(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)
        self.umi_threshold_method = args.umi_threshold_method
        self.min_support_reads = int(args.min_support_reads)
        
        # data
        with open(args.raw_read_count_file) as f:
            self.count_dict = json.load(f)

        self.raw_umi = 0
        self.total_corrected_umi = 0
        self.del_umi = 0
        self.umi_threshold = 1 # if not set explicitly, use 1 as default
        self.umi_count_dict = defaultdict(int)
        self.umi_array = []

        match_dir_dict = utils.parse_match_dir(args.match_dir)
        self.df_tsne = match_dir_dict['df_tsne']

        self.df_filter_tsne = self.df_tsne.copy()

        # out
        self.otsu_plot = f'{self.out_prefix}_otsu.png'
        self.raw_umi_count_file = f'{self.out_prefix}_corrected_UMI_count.json'
        self.filter_tsne_file = f'{self.out_prefix}_filtered_UMI_tsne.csv'
    

    @utils.add_log
    def correct_umi(self):
        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                n_corrected_umi, _n_corrected_read = Count.correct_umi(self.count_dict[barcode][ref])
                if self.debug:
                    print(f'{barcode} {ref} {n_corrected_umi}')
                self.total_corrected_umi += n_corrected_umi
                self.raw_umi += len(self.count_dict[barcode][ref])

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
    def filter_read(self):
        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                for umi in self.count_dict[barcode][ref]:
                    read_count = self.count_dict[barcode][ref][umi]
                    if read_count < self.min_support_reads:
                        self.del_umi += 1
                        self.count_dict[barcode][ref][umi] = 0

        self.add_metric(
            name='Minimum Number of Reads to Support a UMI',
            value=self.min_support_reads,
            help_info='filter UMI with less than this number of reads',
        )

        self.add_metric(
            name='Number of Filtered UMI',
            value=self.del_umi,
            total=self.raw_umi,
            help_info='filter UMI according to `min_support_read`',
        )
        
    
    @utils.add_log
    def set_umi_count_dict(self):
        for barcode in self.count_dict:
            umi_count = 0
            for ref in self.count_dict[barcode]:
                for umi in self.count_dict[barcode][ref]:
                    if self.count_dict[barcode][ref][umi] > 0:
                        umi_count += 1
            self.umi_count_dict[barcode] = umi_count
        if self.debug:
            print(self.umi_count_dict)   

    @utils.add_log
    def set_umi_array(self):
        self.umi_array = list(self.umi_count_dict.values())
                  

    def get_umi_threshold(self):
        if self.umi_threshold_method == 'auto':
            self.auto_threshold()
        elif self.umi_threshold_method == 'otsu':
            self.otsu_threshold()
        elif self.umi_threshold_method == 'hard':
            self.hard_threshold()

        self.add_metric('UMI Threshold Method', self.umi_threshold_method)
        self.add_metric('UMI Threshold', self.umi_threshold)

    @utils.add_log
    def otsu_threshold(self):
        array = np.log2(self.umi_array)
        hist = otsu.array2hist(array)
        thresh = otsu.threshold_otsu(hist)
        otsu.makePlot(hist, thresh, self.otsu_plot)

        self.umi_threshold = int(2 ** thresh)

    @utils.add_log
    def auto_threshold(self):
        """
        threhold = 99 percentile of all cell virus UMIs / 10
        """
        umi_array = self.umi_array

        cell_99th = len(umi_array) // 100
        sorted_umis = sorted(umi_array, reverse=True)
        percentile_99_umi = sorted_umis[cell_99th]
        self.umi_threshold = int(percentile_99_umi / 10)
    
    @utils.add_log
    def hard_threshold(self):
        self.umi_threshold = int(self.args.umi_hard_threshold)

    @utils.add_log
    def filter_umi(self):
        for barcode in self.umi_count_dict:
            if self.umi_count_dict[barcode] < self.umi_threshold:
                self.umi_count_dict[barcode] = 0
    
    @utils.add_log
    def add_umi_write_tsne(self):
        self.df_filter_tsne['UMI'] = pd.Series(self.umi_count_dict)
        self.df_filter_tsne['UMI'].fillna(0, inplace=True)
        self.df_filter_tsne.to_csv(self.filter_tsne_file)


    def add_some_metrics(self):
        df = self.df_filter_tsne

        cell_total = len(df)
        df_virus = df[df['UMI'] > 0]
        cell_with_virus = len(df_virus)
        self.add_metric(
            name='Number of Cells with Virus UMI after Filtering',
            value=cell_with_virus,
            total=cell_total,
        )


    def run(self):
        self.correct_umi()
        self.filter_read()
        self.set_umi_count_dict()
        self.set_umi_array()
        self.get_umi_threshold()
        self.filter_umi()
        self.add_umi_write_tsne()
        self.add_some_metrics()


