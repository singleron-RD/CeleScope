import pandas as pd

from celescope.tools.featureCounts import FeatureCounts as Fc, get_opts_featureCounts as super_opts

class FeatureCounts(Fc):
    def __init__(self, args, display_title=None):
        Fc.__init__(self, args, display_title=display_title)

        self.filter_bam = f'{outdir}/{sample}_filter.bam'

    @staticmethod
    def get_valid_barcode(filter_umi_file):
        df = pd.read_csv


    def filter_bam(self):
        df_barcode
        

        

    def run(self):
        pass



def get_opts_featureCounts(parser, sub_program):
    super_opts(parser, sub_program)
    if sub_program:
        parser.add_argument(
            "--gtf",
            help="Optional. Genome gtf file. Use absolute path or relative path to `genomeDir`.",
        )
        parser.add_argument('--filter_umi_file', help='filter umi file', required=True)
        parser.add_argument('--filter_read_count_json', help='filter_read_count_json file', required=True)