import subprocess
import sys
import pysam

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

    def __init__(self, args, display_title="Filtering"):
        super().__init__(args, display_title)

        self.add_metric(
            "reference allele threshold method",
            self.args.ref_threshold_method,
        )

        self.add_metric(
            "alternate allele threshold method",
            self.args.alt_threshold_method,
        )

        # out
        self.bcftools_vcf_file = f"{self.out_prefix}_bcftools_filtered.vcf"
        self.out_vcf_file = f"{self.out_prefix}_filtered.vcf"

    def bcftools_filter(self):
        cmd = f"bcftools filter -i '{self.args.bcftools_filter}' {self.args.vcf} -o {self.bcftools_vcf_file} -Ov"
        sys.stderr.write(f"{cmd}\n")
        subprocess.check_call(cmd, shell=True)

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
            ad = record.samples[sample]["AD"]
            ref_count_array.append(ad[0])
            alt_count_array.append(ad[1])

        return ref_count_array, alt_count_array

    def get_threshold(self, count_array, threshold_method, otsu_plot_path=None):
        runner = Threshold(
            count_array,
            threshold_method=threshold_method,
            otsu_plot_path=otsu_plot_path,
        )
        threshold = runner.run()
        return threshold

    @staticmethod
    def filter_array(count_array, threshold):
        return [x if x >= threshold else 0 for x in count_array]

    @staticmethod
    def get_genotype(ref_count, alt_count, vaf):
        if ref_count > 0 and alt_count > 0:
            af = alt_count / (alt_count + ref_count)
            if vaf <= af <= 1 - vaf:
                return (0, 1)
            elif af < vaf:
                return (0, 0)
            else:
                return (1, 1)
        elif ref_count > 0 and alt_count == 0:
            return (0, 0)
        elif ref_count == 0 and alt_count > 0:
            return (1, 1)
        elif ref_count == 0 and alt_count == 0:
            return (None, None)

    @utils.add_log
    def run(self):
        self.bcftools_filter()

        if self.args.ref_threshold_method == "otsu":
            ref_otsu_plot_dir = f"{self.outdir}/ref_otsu_plots/"
            utils.check_mkdir(ref_otsu_plot_dir)

        if self.args.alt_threshold_method == "otsu":
            alt_otsu_plot_dir = f"{self.outdir}/alt_otsu_plots/"
            utils.check_mkdir(alt_otsu_plot_dir)

        with pysam.VariantFile(self.bcftools_vcf_file) as vcf_in:
            header = vcf_in.header
            header.add_meta(
                "ref_threshold_method", value=self.args.ref_threshold_method
            )
            header.add_meta(
                "alt_threshold_method", value=self.args.alt_threshold_method
            )
            header.add_meta(
                "INFO",
                items=[
                    ("ID", "REF_T"),
                    ("Number", 1),
                    ("Type", "Integer"),
                    ("Description", "Reference allele count threshold"),
                ],
            )
            header.add_meta(
                "INFO",
                items=[
                    ("ID", "ALT_T"),
                    ("Number", 1),
                    ("Type", "Integer"),
                    ("Description", "Alternate allele count threshold"),
                ],
            )
            with pysam.VariantFile(
                self.out_vcf_file, "w", header=vcf_in.header
            ) as vcf_out:
                for record in vcf_in.fetch():
                    name = f"{record.chrom}_{record.pos}"
                    ref_count_array, alt_count_array = self.get_count_array(record)
                    ref_otsu_plot_path = (
                        f"{ref_otsu_plot_dir}/{name}_ref.pdf"
                        if self.args.ref_threshold_method == "otsu"
                        else None
                    )
                    ref_threshold = self.get_threshold(
                        ref_count_array,
                        self.args.ref_threshold_method,
                        ref_otsu_plot_path,
                    )
                    alt_otsu_plot_path = (
                        f"{alt_otsu_plot_dir}/{name}_alt.pdf"
                        if self.args.alt_threshold_method == "otsu"
                        else None
                    )
                    alt_threshold = self.get_threshold(
                        alt_count_array,
                        self.args.alt_threshold_method,
                        alt_otsu_plot_path,
                    )
                    ref_threshold = max(ref_threshold, self.args.ref_min_support_read)
                    alt_threshold = max(alt_threshold, self.args.alt_min_support_read)
                    ref_filtered_count_array = self.filter_array(
                        ref_count_array, ref_threshold
                    )
                    alt_filtered_count_array = self.filter_array(
                        alt_count_array, alt_threshold
                    )

                    new_record = record.copy()
                    new_record.info.__setitem__("REF_T", ref_threshold)
                    new_record.info.__setitem__("ALT_T", alt_threshold)
                    for index, sample in enumerate(record.samples):
                        ref_count = ref_filtered_count_array[index]
                        alt_count = alt_filtered_count_array[index]
                        new_record.samples[sample]["AD"] = (ref_count, alt_count)
                        genotype = self.get_genotype(
                            ref_count, alt_count, self.args.vaf
                        )
                        new_record.samples[sample]["GT"] = genotype

                    vcf_out.write(new_record)


def filter_snp(args):
    with Filter_snp(args) as runner:
        runner.run()


def get_opts_filter_snp(parser, sub_program):
    parser.add_argument(
        "--bcftools_filter",
        help="bcftools filter -i expression",
        default="QUAL>=100",
    )
    parser.add_argument(
        "--ref_threshold_method",
        default="none",
        choices=["otsu", "none"],
        help=HELP_DICT["threshold_method"],
    )
    parser.add_argument(
        "--alt_threshold_method",
        default="none",
        choices=["otsu", "none"],
        help=HELP_DICT["threshold_method"],
    )
    parser.add_argument(
        "--ref_min_support_read",
        type=int,
        help="minimum supporting read number for ref.",
        default=2,
    )
    parser.add_argument(
        "--alt_min_support_read",
        type=int,
        help="minimum supporting read number for alt.",
        default=2,
    )
    parser.add_argument(
        "--vaf",
        help="VAF threshold. For each cell, if VAF_threshold <= ALT/DP <= 1 - VAF_threshold, then the genotype is called heterozygous; otherwise, it is called homozygous.",
        type=float,
        default=0.2,
    )
    if sub_program:
        parser.add_argument("--vcf", help="norm vcf file")
        s_common(parser)
