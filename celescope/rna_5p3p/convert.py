"""
convert 5-prime and 3-prime to the same pattern.
"""

from celescope.tools import utils, parse_chemistry
from celescope.tools.step import Step, s_common
import subprocess
import pysam


class Convert(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fq1_5p = args.fq1_5p.split(",")
        self.fq1_3p = args.fq1_3p.split(",")
        self.fq2_5p = args.fq2_5p.split(",")
        self.fq2_3p = args.fq2_3p.split(",")

        chemistry_version = args.chemistry.split("-")[1]
        chemistry_5p, chemistry_3p = (
            f"rna_5p-{chemistry_version}",
            f"rna_3p-{chemistry_version}",
        )
        chemistry_dict = parse_chemistry.get_chemistry_dict()
        self.chemistry_dict_5p, self.chemistry_dict_3p = (
            chemistry_dict[chemistry_5p],
            chemistry_dict[chemistry_3p],
        )

    def write_3p(self):
        for i, (fn1, fn2) in enumerate(zip(self.fq1_3p, self.fq2_3p), start=1):
            out_fn1 = f"{self.out_prefix}_3p{i}_R1.fq.gz"
            out_fn2 = f"{self.out_prefix}_3p{i}_R2.fq.gz"
            fh_3p = utils.generic_open(out_fn1, "wt")
            with pysam.FastxFile(fn1, persist=False) as fq:
                for entry1 in fq:
                    name1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    bc_list, bc_quality_list, umi, umi_quality = (
                        parse_chemistry.get_raw_umi_bc_and_quality(
                            seq1, qual1, self.chemistry_dict_3p["pattern_dict"]
                        )
                    )
                    bc = "".join(bc_list)
                    bc_quality = "".join(bc_quality_list)
                    fh_3p.write(
                        utils.fastq_line(name1, bc + umi, bc_quality + umi_quality)
                    )
            cmd = f"ln -s -f {fn2} {out_fn2}"
            subprocess.check_call(cmd, shell=True)

    def write_5p(self):
        for i, (fn1, fn2) in enumerate(zip(self.fq1_5p, self.fq2_5p), start=1):
            out_fn1 = f"{self.out_prefix}_5p{i}_R1.fq.gz"
            out_fn2 = f"{self.out_prefix}_5p{i}_R2.fq.gz"
            fh_5p_1 = utils.generic_open(out_fn1, "wt")
            fh_5p_2 = utils.generic_open(out_fn2, "wt")
            with pysam.FastxFile(fn1, persist=False) as fq:
                for entry1 in fq:
                    name1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    bc_list, bc_quality_list, umi, umi_quality = (
                        parse_chemistry.get_raw_umi_bc_and_quality(
                            seq1,
                            qual1,
                            self.chemistry_dict_5p["pattern_dict"],
                            reverse_complement=True,
                        )
                    )
                    bc = "".join(bc_list)
                    bc_quality = "".join(bc_quality_list)
                    fh_5p_1.write(
                        utils.fastq_line(name1, bc + umi, bc_quality + umi_quality)
                    )
            fh_5p_1.close()
            with pysam.FastxFile(fn2, persist=False) as fq:
                for entry2 in fq:
                    name2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    seq2 = utils.reverse_complement(seq2)
                    qual2 = qual2[::-1]
                    fh_5p_2.write(utils.fastq_line(name2, seq2, qual2))
            fh_5p_2.close()

    def run(self):
        self.write_3p()
        self.write_5p()


@utils.add_log
def convert(args):
    with Convert(args) as runner:
        runner.run()


def get_opts_convert(parser, sub_program=True):
    parser.add_argument(
        "--chemistry",
        required=True,
        choices=["5p3p-1", "5p3p-2", "5p3p-3"],
        help="rna 5p3p chemistry version",
    )
    if sub_program:
        parser.add_argument(
            "--fq1_3p",
            help="3 prime R1 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser.add_argument(
            "--fq2_3p",
            help="3 prime R2 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser.add_argument(
            "--fq1_5p",
            help="5 prime R1 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser.add_argument(
            "--fq2_5p",
            help="5 prime R2 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser = s_common(parser)

    return parser
