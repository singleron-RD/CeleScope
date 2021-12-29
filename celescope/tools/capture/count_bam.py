import json

import pysam
import numpy as np

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT, HELP_INFO_DICT


def get_opts_count_bam(parser, sub_program):
    parser.add_argument("--min_query_length", help='Minimum query length.', default=35)

    if sub_program:
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--capture_bam', required=True)
        s_common(parser)


class Count_bam(Step):

    def __init__(self, args, display_title='Count'):
        super().__init__(args, display_title)

        # set
        self.min_query_length = int(args.min_query_length)
        self.capture_bam = args.capture_bam

        # read barcodes
        match_dir_dict = utils.parse_match_dir(args.match_dir)
        self.match_barcode = match_dir_dict['match_barcode']
        self.n_match_barcode = match_dir_dict['n_match_barcode']
        self.add_metric(
            name=HELP_INFO_DICT['matched_barcode_number']['display'],
            value=self.n_match_barcode,
            help_info=HELP_INFO_DICT['matched_barcode_number']['info']
        )

        # data
        self.total_corrected_umi = 0
        self.count_dict = utils.genDict(dim=3)

        # out
        self.raw_read_count_file = f'{self.out_prefix}_raw_read_count.json'

    @utils.add_log
    def process_bam(self):
        # process bam
        samfile = pysam.AlignmentFile(self.args.capture_bam, "rb")
        for read in samfile:
            ref = read.reference_name
            query_length = read.infer_query_length()
            attr = read.query_name.split('_')
            barcode = attr[0]
            umi = attr[1]
            if (barcode in self.match_barcode) and (query_length >= self.min_query_length):
                self.count_dict[barcode][ref][umi] += 1

    @utils.add_log
    def add_some_metrics(self):
        read_count_list = []
        for barcode in self.count_dict:
            for ref in self.count_dict[barcode]:
                for umi in self.count_dict[barcode][ref]:
                    read_count_list.append(self.count_dict[barcode][ref][umi])
        mean_read_count_per_umi = round(float(np.mean(read_count_list)), 2)
        self.add_metric(
            name='Mead Read Count per UMI',
            value=mean_read_count_per_umi,
            help_info='you can use this value to determine `min_support_read`'
        )

        n_positive_cell = len(self.count_dict)
        self.add_metric(
            name='Number of positive cells',
            value=n_positive_cell,
            total=self.n_match_barcode,
            help_info='this number is calculated on raw read counts(i.e. before filtering)'
        )

    @utils.add_log
    def write_count_file(self):

        with open(self.raw_read_count_file, 'w') as fp:
            json.dump(self.count_dict, fp, indent=4)

    @utils.add_log
    def run(self):
        self.process_bam()
        self.add_some_metrics()
        self.write_count_file()


