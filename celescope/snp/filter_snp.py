import pysam
import subprocess

from celescope.tools import utils
from celescope.tools.capture.threshold import Threshold
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT

class Filter_snp(Step):
    """
    ## Features
    - Filter out `ref` and `alt` alleles that do not have enough reads to support.

    ## Output
    - `{sample}_test1_filtered.vcf` VCF file after filtering. Alleles read counts that do not have enough reads to support are set to zero. 
    Genotypes are changed accordingly.
    """
    def __init__(self, args, display_title='Filtering'):
        super().__init__(args, display_title)
        self.vcf = args.vcf
        self.threshold_method = args.threshold_method
        self.hard_threshold = args.hard_threshold

        self.add_metric(
            'threshold method',
            self.threshold_method,
        )

        # out
        self.out_vcf_file = f'{self.out_prefix}_filtered.vcf'


    @staticmethod
    def get_count_array(record):
        """
        get allele count array

        Args
            record: pysam vcf record
            entry: str, 'ref' or 'alt'
        Return
            (ref_count_array, alt_count_array)
        """
        ref_count_array = []
        alt_count_array = []
        for sample in record.samples:
            ad = record.samples[sample]['AD']
            ref_count_array.append(ad[0])
            alt_count_array.append(ad[1])

        return ref_count_array, alt_count_array

    def get_threshold(self, count_array):
        runner = Threshold(
            count_array,
            threshold_method=self.threshold_method,
            hard_threshold=self.hard_threshold,
        )
        threshold = runner.run()
        return threshold

    @staticmethod
    def filter_array(count_array, threshold):
        return [x if x >= threshold else 0 for x in count_array]

    @staticmethod
    def get_genotype(ref_count, alt_count):
        if ref_count > 0 and alt_count > 0:
            return (0,1)
        elif ref_count > 0 and alt_count == 0:
            return (0,0)
        elif ref_count == 0 and alt_count > 0:
            return (1,1)
        elif ref_count == 0 and alt_count == 0:
            return (None, None)

    @utils.add_log
    def run(self):
        if self.threshold_method == 'none':

            cmd = f'cp {self.vcf} {self.out_vcf_file}'
            subprocess.check_call(cmd, shell=True)
            return

        with pysam.VariantFile(self.vcf) as vcf_in:
            header = vcf_in.header
            header.add_meta('threshold_method', value=self.threshold_method)
            header.add_meta('INFO', items=[('ID',"REF_T"), ('Number',1), ('Type','Integer'), ('Description','Reference allele count threshold')])
            header.add_meta('INFO', items=[('ID',"ALT_T"), ('Number',1), ('Type','Integer'), ('Description','Alternate allele count threshold')])
            with pysam.VariantFile(self.out_vcf_file, 'w', header=vcf_in.header) as vcf_out:
                for record in vcf_in.fetch():
                    ref_count_array, alt_count_array = self.get_count_array(record)
                    ref_threshold = self.get_threshold(ref_count_array)
                    alt_threshold = self.get_threshold(alt_count_array)
                    ref_filtered_count_array = self.filter_array(ref_count_array, ref_threshold)
                    alt_filtered_count_array = self.filter_array(alt_count_array, alt_threshold)

                    new_record = record.copy()
                    new_record.info.__setitem__('REF_T', ref_threshold)
                    new_record.info.__setitem__('ALT_T', alt_threshold)
                    for index, sample in enumerate(record.samples):
                        ref_count = ref_filtered_count_array[index]
                        alt_count = alt_filtered_count_array[index]
                        new_record.samples[sample]['AD'] = (ref_count, alt_count)
                        genotype = self.get_genotype(ref_count, alt_count)
                        new_record.samples[sample]['GT'] = genotype

                    vcf_out.write(new_record)


def filter_snp(args):
    with Filter_snp(args) as runner:
        runner.run()


def get_opts_filter_snp(parser, sub_program):
    parser.add_argument('--threshold_method', default='auto', choices=['otsu', 'auto', 'hard', 'none'], help=HELP_DICT['threshold_method'])
    parser.add_argument(
        "--hard_threshold",
        help='int, use together with `--threshold_method hard`',
    )
    if sub_program:
        parser.add_argument("--vcf", help="norm vcf file")
        s_common(parser)




        





