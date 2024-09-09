import json
import collections

import pandas as pd
import pysam

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT


class FeatureCounts(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.filter_bam = f"{self.outdir}/{self.sample}_filter.bam"
        self.coord_sorted_bam = (
            f"{self.outdir}/{self.sample}_filter.bam.featureCounts.bam"
        )
        self.name_sorted_bam = f"{self.outdir}/{self.sample}_filter_name_sorted.bam"
        self.log_file = f"{self.outdir}/{self.sample}.summary"

    @staticmethod
    def get_valid_barcodes(filter_umi_file):
        df = pd.read_csv(filter_umi_file)
        df = df[df["sum_UMI"] > 0]
        valid_barcodes = set(df["barcode"])
        return valid_barcodes

    @staticmethod
    def get_valid_umis(filter_read_count_json, valid_barcodes):
        barcode_umis = collections.defaultdict(set)
        with open(filter_read_count_json) as fh:
            dic = json.load(fh)
        for barcode in dic:
            if barcode in valid_barcodes:
                for ref in dic[barcode]:
                    for umi in dic[barcode][ref]:
                        if dic[barcode][ref][umi] > 0:
                            barcode_umis[barcode].add(umi)
        return barcode_umis

    def run_filter(self):
        valid_barcodes = FeatureCounts.get_valid_barcodes(self.args.filter_umi_file)
        barcode_umis = FeatureCounts.get_valid_umis(
            self.args.filter_read_count_json, valid_barcodes
        )

        with pysam.AlignmentFile(self.args.bam, "rb") as raw_bam:
            header = raw_bam.header
            with pysam.AlignmentFile(self.filter_bam, "w", header=header) as filter_bam:
                for read in raw_bam:
                    attr = read.query_name.split(":")
                    barcode = attr[0]
                    umi = attr[1]
                    if barcode in valid_barcodes and umi in barcode_umis[barcode]:
                        filter_bam.write(read)

    @utils.add_log
    def run_featureCounts(self):
        cmd = (
            "featureCounts "
            "-s 1 "
            f"-a {self.args.gtf} "
            f"-o {self.outdir}/{self.sample} "  # not bam
            "-R BAM "
            f"-T {self.thread} "
            f"-t gene "
            f"{self.filter_bam} "
            "2>&1 "
        )
        if self.args.featureCounts_param:
            cmd += " " + self.args.featureCounts_param
        self.debug_subprocess_call(cmd)

        utils.sort_bam(
            input_bam=self.coord_sorted_bam, output_bam=self.name_sorted_bam, by="name"
        )

    @staticmethod
    def read_log(log_file):
        """
        Args:
            log_file: featureCounts log summary file
        Returns:
            log_dict: {'Assigned': 123, ...}
        """
        # skip first line
        df = pd.read_csv(
            log_file, sep="\t", header=None, names=["name", "value"], skiprows=1
        )
        log_dict = df.set_index("name")["value"].to_dict()
        return log_dict

    def add_metrics(self):
        log_dict = self.read_log(self.log_file)
        total = sum(log_dict.values())
        self.add_metric(
            name="Reads Assigned",
            value=log_dict["Assigned"],
            total=total,
            help_info="Reads that can be successfully assigned",
        )

        self.add_metric(
            name="Reads Unassigned_NoFeatures",
            value=log_dict["Unassigned_NoFeatures"],
            total=total,
            help_info="Reads that can not be ssigned",
        )

        self.add_metric(
            name="Reads Unassigned_Ambiguity",
            value=log_dict["Unassigned_Ambiguity"],
            total=total,
            help_info="Reads that overlap two or more features",
        )

    def run(self):
        self.run_filter()
        self.run_featureCounts()
        self.add_metrics()


@utils.add_log
def featureCounts(args):
    if args.gtf:
        with FeatureCounts(args) as runner:
            runner.run()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument(
        "--gtf",
        help="Optional. Genome gtf file. Use absolute path or relative path to `genomeDir`.",
    )
    if sub_program:
        parser.add_argument("--bam", help="input bam file", required=True)
        parser.add_argument("--filter_umi_file", help="filter umi file", required=True)
        parser.add_argument(
            "--filter_read_count_json",
            help="filter_read_count_json file",
            required=True,
        )
        parser.add_argument(
            "--featureCounts_param", help=HELP_DICT["additional_param"], default=""
        )
        s_common(parser)
