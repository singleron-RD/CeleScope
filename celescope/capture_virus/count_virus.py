import pandas as pd
import pysam
import numpy as np

import celescope.tools.utils as utils
from celescope.tools.step import Step,s_common
from celescope.tools.count import Count
from celescope.__init__ import HELP_INFO_DICT


def get_opts_count_virus(parser, sub_program):
    parser.add_argument("--min_query_length", help='Minimum query length.', default=35)

    if sub_program:
        parser.add_argument('--match_dir', help='matched rna_virus directory', required=True)
        parser.add_argument('--virus_bam', required=True)
        s_common(parser)

class Count_virus(Step):

    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        # set
        self.min_query_length = int(args.min_query_length)

        # read barcodes
        match_dir_dict = utils.parse_match_dir(args.match_dir)
        self.match_barcode = match_dir_dict['match_barcode']
        self.n_match_barcode = match_dir_dict['n_match_barcode']
        self.df_tsne = match_dir_dict['df_tsne']
        self.add_metric(
            name=HELP_INFO_DICT['matched_barcode_number']['display'],
            value=self.n_match_barcode,
            help_info=HELP_INFO_DICT['matched_barcode_number']['info']
        )

        # data
        self.total_corrected_umi = 0
        self.count_dic = utils.genDict(dim=3)
        self.df_umi = pd.DataFrame()

        # out
        self.raw_read_count_file = f'{self.out_prefix}_raw_read_count.csv'
        self.raw_umi_count_file = f'{self.out_prefix}_raw_UMI_count.csv'
        self.umi_tsne_file = f'{self.out_prefix}_UMI_tsne.csv'


    @utils.add_log
    def process_bam(self):
        # process bam
        samfile = pysam.AlignmentFile(self.args.virus_bam, "rb")
        for read in samfile:
            ref = read.reference_name
            query_length = read.infer_query_length()
            attr = read.query_name.split('_')
            barcode = attr[0]
            umi = attr[1]
            if (barcode in self.match_barcode) and (query_length >= self.min_query_length):
                self.count_dic[barcode][ref][umi] += 1

    @utils.add_log
    def correct_umi(self):
        for barcode in self.count_dic:
            for ref in self.count_dic[barcode]:
                n_corrected_umi, _n_corrected_read = Count.correct_umi(self.count_dic[barcode][ref])
                self.total_corrected_umi += n_corrected_umi
        self.add_metric(
            name='Corrected UMI Number',
            value=self.total_corrected_umi,
            help_info='number of corrected UMI',
        )
                

    @utils.add_log
    def write_count_file(self):

        # write dic to pandas df
        rows = []
        for barcode in self.count_dic:
            for ref in self.count_dic[barcode]:
                for umi in self.count_dic[barcode][ref]:
                    read_count = self.count_dic[barcode][ref][umi]
                    rows.append([barcode, ref, umi, read_count])

        df_read = pd.DataFrame(
            rows,
            columns=[
                "barcode",
                "ref",
                "UMI",
                "read_count"]
        )
        df_read.to_csv(self.raw_read_count_file, index=False)

        self.df_umi = df_read.groupby(["barcode", "ref"]).agg({'UMI': 'count'}).reset_index()
        self.df_umi.to_csv(self.raw_umi_count_file)

    
    @utils.add_log
    def write_tsne(self):
        '''
        aggregate UMI count and write to csv
        '''
        df_sum_umi = self.df_umi.groupby(["barcode"]).agg({'UMI': 'sum'}).reset_index()
        df_umi_tsne = pd.merge(self.df_tsne, df_sum_umi, on="barcode", how="left")
        df_umi_tsne.index = df_umi_tsne.barcode
        df_umi_tsne.drop(columns=["barcode"], inplace=True)
        df_umi_tsne.to_csv(self.umi_tsne_file)

    def run(self):
        self.process_bam()
        self.correct_umi()
        self.write_count_file()
        self.write_tsne()


@utils.add_log
def count_virus(args):

    with Count_virus(args, display_title='Count') as runner:
        runner.run()




