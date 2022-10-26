import collections
import subprocess
import pandas as pd
import pysam

from celescope.tools.featureCounts import FeatureCounts as Fc
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT


class FeatureCounts_capture_virus(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        
        self.filter_bam = f'{self.outdir}/{self.sample}_filter.bam'
        self.coord_sorted_bam = f'{self.outdir}/{self.sample}_filter.bam.featureCounts.bam'
        self.name_sorted_bam = f'{self.outdir}/{self.sample}_filter_name_sorted.bam'
        self.log_file = f'{self.outdir}/{self.sample}.summary'
    
    @staticmethod
    def get_valid_barcodes(filter_umi_file):
        df = pd.read_table(filter_umi_file)
        df = df[df['UMI'] > 0]
        valid_barcodes = set(df['barcode'])
        return valid_barcodes    
    
    @staticmethod
    def get_valid_umis(filter_read_count_df, valid_barcodes):
        df = pd.read_table(filter_read_count_df)
        filter_bc_df = df[df["barcode"].isin(valid_barcodes)][["barcode","UMI"]].drop_duplicates()
        barcode_umis = collections.defaultdict(set)
        all_bc = filter_bc_df["barcode"].to_list()
        for bc in all_bc:
            umi_set = set(filter_bc_df[filter_bc_df["barcode"]==bc]["UMI"].to_list())
            barcode_umis[bc] = barcode_umis[bc] | umi_set
        return barcode_umis

    def run_filter(self):
        valid_barcodes = FeatureCounts_capture_virus.get_valid_barcodes(self.args.filter_umi_file)
        barcode_umis = FeatureCounts_capture_virus.get_valid_umis(self.args.filter_read_count_json, valid_barcodes)

        with pysam.AlignmentFile(self.args.bam, "rb") as raw_bam:
            header = raw_bam.header
            with pysam.AlignmentFile(self.filter_bam, "w", header=header) as filter_bam:
                for read in raw_bam:
                    attr = read.query_name.split('_')
                    barcode = attr[0]
                    umi = attr[1]
                    if barcode in valid_barcodes and umi in barcode_umis[barcode]:
                        filter_bam.write(read)      

    @utils.add_log
    def run_featureCounts(self):
        cmd = (
            'featureCounts '
            '-s 1 '
            f'-a {self.args.gtf} '
            f'-o {self.outdir}/{self.sample} '  # not bam
            '-R BAM '
            f'-T {self.thread} '
            f'-t gene '
            f'{self.filter_bam} '
            '2>&1 '
        )
        if self.args.featureCounts_param:
            cmd += (" " + self.args.featureCounts_param)
        subprocess.check_call(cmd, shell=True)
        
        utils.sort_bam(
            input_bam = self.coord_sorted_bam,
            output_bam = self.name_sorted_bam,
            by='name'
        )

    def add_metrics(self):
        log_dict = Fc.read_log(self.log_file)
        total = sum(log_dict .values())
        self.add_metric(
            name='Reads Assigned',
            value=log_dict['Assigned'],
            total=total
        )

        self.add_metric(
            name='Reads Unassigned_NoFeatures',
            value=log_dict['Unassigned_NoFeatures'],
            total=total
        )

        self.add_metric(
            name='Reads Unassigned_Ambiguity',
            value=log_dict['Unassigned_Ambiguity'],
            total=total
        )

    def run(self):
        self.run_filter()
        self.run_featureCounts()
        self.add_metrics()

@utils.add_log
def featureCounts_capture_virus(args):

    step_name = 'featureCounts_capture_virus'
    runner = FeatureCounts_capture_virus(args, step_name)
    runner.run()


def get_opts_featureCounts_capture_virus(parser, sub_program):
    parser.add_argument(
        "--gtf",
        help="Optional. Genome gtf file. Use absolute path or relative path to `genomeDir`.",
        )
    if sub_program:
        parser.add_argument('--bam', help='input bam file', required=True)
        parser.add_argument('--filter_umi_file', help='filter umi file', required=True)
        parser.add_argument('--filter_read_count_json', help='filter_read_count_json file', required=True)
        parser.add_argument('--featureCounts_param', help=HELP_DICT['additional_param'], default="")
        s_common(parser)