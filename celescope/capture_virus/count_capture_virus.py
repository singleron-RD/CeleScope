import pandas as pd
import pysam

import celescope.tools.utils as utils
from celescope.tools.step import Step,s_common

class Count_capture_virus(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # set
        self.min_query_length = int(args.min_query_length)
        self.min_support_read = int(args.min_support_read)

        # read barcodes
        self.match_cell_barcodes, _match_cell_number = utils.read_barcode_file(args.match_dir)

        # out
        self.out_read_count_file = f'{self.out_prefix}_virus_read_count.tsv'
        self.out_umi_count_file = f'{self.out_prefix}_virus_UMI_count.tsv'

    @utils.add_log
    def sum_virus(self):
        # process bam
        samfile = pysam.AlignmentFile(self.args.virus_bam, "rb")
        count_dic = utils.genDict(dim=3)
        for read in samfile:
            tag = read.reference_name
            query_length = read.infer_query_length()
            attr = read.query_name.split('_')
            barcode = attr[0]
            umi = attr[1]
            if (barcode in self.match_cell_barcodes) and (query_length >= self.min_query_length):
                count_dic[barcode][tag][umi] += 1

        # write dic to pandas df
        rows = []
        for barcode in count_dic:
            for tag in count_dic[barcode]:
                for umi in count_dic[barcode][tag]:
                    rows.append([barcode, tag, umi, count_dic[barcode][tag][umi]])

        df_read = df_read = pd.DataFrame(
            rows,
            columns=[
                "barcode",
                "tag",
                "UMI",
                "read_count"])
        df_read.to_csv(self.out_read_count_file, sep="\t", index=False)

        df_valid = df_read[df_read["read_count"] >= self.min_support_read].groupby(
            ["barcode", "tag"]).agg({"UMI": "count"})
        if df_valid.shape[0] == 0:
            self.sum_virus.logger.warning("No cell virus UMI found!")

        df_umi = df_valid.groupby(["barcode", "tag"]).agg({"UMI": "count"})
        df_umi.to_csv(self.out_umi_count_file, sep="\t")

    def run(self):
        self.sum_virus()


@utils.add_log
def count_capture_virus(args):

    step_name = 'count_capture_virus'
    runner = Count_capture_virus(args, step_name)
    runner.run()


def get_opts_count_capture_virus(parser, sub_program):
    parser.add_argument("--min_query_length", help='Minimum query length.', default=35)
    parser.add_argument("--min_support_read", help='Minimum number of reads supporting a UMI', default=2)
    if sub_program:
        parser.add_argument('--match_dir', help='matched rna_virus directory', required=True)
        parser.add_argument('--virus_bam', required=True)
        s_common(parser)
