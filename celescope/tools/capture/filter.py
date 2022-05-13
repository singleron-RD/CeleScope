import json

import pandas as pd

from celescope.tools import utils
from celescope.tools.count import Count
from celescope.tools.step import Step, s_common
from celescope.tools.capture.threshold import Threshold
from celescope.__init__ import HELP_DICT
from celescope.tools.capture.__init__ import SUM_UMI_COLNAME



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
        default='otsu'
    )
    parser.add_argument(
        "--umi_hard_threshold",
        help='int, use together with `--umi_threshold_method hard`',
    )
    parser.add_argument(
        "--auto_coef",
        help='int, threshold = top 1% positive cell count / auto_coef',
        default=3,
    )
    parser.add_argument(
        "--otsu_log_base",
        help='raw counts are first log transformed before thresholding. This argument is the log base. Commonly used values are 2 and 10.',
        default=10,
        type=int,
    )


    if sub_program:
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--raw_read_count_file', required=True)
        s_common(parser)


class Filter(Step):
    """    
    ## Features
    - Correct single-base errors in UMIs due to sequencing, amplification, etc.
    - Filter background UMIs base on a UMI threshold.
    There are three methods to determine the UMI threshold:
        - 'auto' : Using a method similar to cell calling method.
        - 'otsu' : UMI counts are first log transformed and then the threshold is determined by [Otsu's method](https://en.wikipedia.org/wiki/Otsu%27s_method)
        - 'hard' : Using User provided UMI threshold.

    ## Output
    - `{sample}_corrected_read_count.json` Read counts after UMI correction.
    - `{sample}_filtered_read_count.json` Filtered read counts.
    - `{sample}_filtered_UMI.csv` Filtered UMI counts.

    """
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

        match_dir_dict = utils.parse_match_dir(args.match_dir)
        self.match_barcode = match_dir_dict['match_barcode']
        self.df_filter_umi = pd.DataFrame(index=list(self.match_barcode)).rename_axis('barcode')

        # out
        self.corrected_read_count_file = f'{self.out_prefix}_corrected_read_count.json'
        self.filter_read_count_file = f'{self.out_prefix}_filtered_read_count.json'
        self.filter_umi_file = f'{self.out_prefix}_filtered_UMI.csv'

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
        self.add_metric(
            'Read Threshold Method', 
            self.args.read_threshold_method,
            help_info=f"""There are three methods to determine the threshold:<br>
1. 'auto' : Using a method similar to cell calling method. threshold = top 1% positive cell count / auto_coef. auto_coef = {self.args.auto_coef} <br>
2. 'otsu' : UMI counts are first log transformed and then the threshold is determined by Otsu's method. otsu_log_base = {self.args.otsu_log_base} <br>
3. 'hard' : Using User provided UMI threshold."""
        )

        if self.args.read_threshold_method == 'auto':
            self.add_metric(
                name='Read auto coeffient',
                value=self.args.auto_coef,
                help_info='threshold = top 1% positive cell count / auto_coef',
            )

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
                hard_threshold=self.args.read_hard_threshold,
                coef=self.args.auto_coef,
                log_base=self.args.otsu_log_base,
            )
            read_threshold = runner.run()
            
            self.read_threshold_dict[ref] = read_threshold
            self.add_metric(
                f'{ref} Read Threshold', 
                read_threshold,
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
        self.add_metric(
            'UMI Threshold Method', 
            self.args.umi_threshold_method,
            help_info="""Use the same threshold method as read count threshold."""
        )

        if self.args.umi_threshold_method == 'auto':
            self.add_metric(
            name='UMI auto coeffient',
            value=self.args.auto_coef,
            help_info='threshold = top 1% positive cell count / auto_coef',
        )

        for ref in self.ref_barcode_umi_dict:
            umi_array = list(self.ref_barcode_umi_dict[ref].values())
            otsu_plot_path = f'{self.out_prefix}_{ref}_UMI_otsu.png'
            runner = Threshold(
                umi_array, 
                threshold_method=self.args.umi_threshold_method, 
                otsu_plot_path=otsu_plot_path,
                hard_threshold=self.args.umi_hard_threshold,
                coef=self.args.auto_coef,
                log_base=self.args.otsu_log_base,
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
    def add_umi_write_csv(self):
        for ref in self.umi_threshold_dict:
            self.df_filter_umi[ref] = pd.Series(self.ref_barcode_umi_dict[ref])
            self.df_filter_umi[ref].fillna(0, inplace=True)

        refs = list(self.umi_threshold_dict.keys()) 
        self.df_filter_umi[SUM_UMI_COLNAME] = self.df_filter_umi[refs].sum(axis=1)
        self.df_filter_umi.to_csv(self.filter_umi_file)


    def add_some_metrics(self):
        df = self.df_filter_umi

        cell_total = len(df)
        df_positive = df[df[SUM_UMI_COLNAME] > 0]
        n_cell_positive = len(df_positive)
        self.add_metric(
            name='Number of Positive Cells after Filtering',
            value=n_cell_positive,
            total=cell_total,
            help_info='Cells with sum_UMI > 0 after filtering',
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

        self.add_umi_write_csv()
        self.add_some_metrics()
