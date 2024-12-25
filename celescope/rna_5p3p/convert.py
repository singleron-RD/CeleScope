"""
1. 读入5p R1,转成3p barcode + UMI
2. 读入3p R1,生成3p barcode + UMI
"""

import subprocess

from xopen import xopen
import pysam

from celescope.tools import utils
from celescope.tools.__init__ import PATTERN_DICT
from celescope.tools.step import Step, s_common
from celescope.tools.barcode import Chemistry, Barcode as Bc


class Convert(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.fq1_5p = args.fq1_5p.split(",")
        self.fq1_3p = args.fq1_3p.split(",")
        self.fq2_5p = args.fq2_5p.split(",")
        self.fq2_3p = args.fq2_3p.split(",")
        if args.chemistry == "5p3p-1":
            self.pattern_dict_5p = Bc.parse_pattern(PATTERN_DICT["rna_5p"])
            self.pattern_dict_3p = Bc.parse_pattern(PATTERN_DICT["rna_3p"])
        elif args.chemistry == "5p3p-2":
            self.pattern_dict_5p = Bc.parse_pattern(PATTERN_DICT["rna_5p.1"])
            self.pattern_dict_3p = Bc.parse_pattern(PATTERN_DICT["rna_3p.1"])
        elif args.chemistry == "5p3p-3":
            self.pattern_dict_5p = Bc.parse_pattern(PATTERN_DICT["rna_5p.2"])
            self.pattern_dict_3p = Bc.parse_pattern(PATTERN_DICT["rna_3p.2"])
        else:
            raise ValueError(
                f"Invalid chemistry {args.chemistry}. chemistry must be one of 5p3p-1 or 5p3p-2"
            )
        whitelist_files = Chemistry.get_whitelist(args.chemistry)
        self.barcode_set_list, self.barcode_mismatch_list = Bc.parse_whitelist_file(
            whitelist_files, n_pattern=len(self.pattern_dict_5p["C"]), n_mismatch=1
        )

    def convert_5p_R1(self, seq1):
        # convert 5p R1 to 3p bc
        seq_list = Bc.get_seq_list(seq1, self.pattern_dict_5p, "C")
        seq_list = [utils.reverse_complement(seq) for seq in seq_list[::-1]]
        bool_valid, bool_corrected, corrected_seq_list = Bc.check_seq_mismatch(
            seq_list, self.barcode_set_list, self.barcode_mismatch_list
        )
        umi = Bc.get_seq_str(seq1, self.pattern_dict_5p["U"])
        umi = utils.reverse_complement(umi)
        if bool_valid:
            bc = "".join(corrected_seq_list)
        else:
            bc = "".join(seq_list)
        return bc, umi

    def write_3p(self):
        for i, (fn1, fn2) in enumerate(zip(self.fq1_3p, self.fq2_3p), start=1):
            out_fn1 = f"{self.out_prefix}_3p{i}_R1.fq.gz"
            out_fn2 = f"{self.out_prefix}_3p{i}_R2.fq.gz"
            fh_3p = xopen(out_fn1, "w")
            with pysam.FastxFile(fn1, persist=False) as fq:
                for entry1 in fq:
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    bc = Bc.get_seq_str(seq1, self.pattern_dict_3p["C"])
                    umi = Bc.get_seq_str(seq1, self.pattern_dict_3p["U"])
                    bc_qual = Bc.get_seq_str(qual1, self.pattern_dict_3p["C"])
                    umi_qual = Bc.get_seq_str(qual1, self.pattern_dict_3p["U"])
                    fh_3p.write(
                        "@{}\n{}\n+\n{}\n".format(header1, bc + umi, bc_qual + umi_qual)
                    )
            cmd = f"ln -s -f {fn2} {out_fn2}"
            subprocess.check_call(cmd, shell=True)

    def write_5p(self):
        for i, (fn1, fn2) in enumerate(zip(self.fq1_5p, self.fq2_5p), start=1):
            out_fn1 = f"{self.out_prefix}_5p{i}_R1.fq.gz"
            out_fn2 = f"{self.out_prefix}_5p{i}_R2.fq.gz"
            fh_5p_1 = xopen(out_fn1, "w")
            fh_5p_2 = xopen(out_fn2, "w")
            with pysam.FastxFile(fn1, persist=False) as fq:
                for entry1 in fq:
                    header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    bc, umi = self.convert_5p_R1(seq1)
                    bc_qual = Bc.get_seq_str(qual1, self.pattern_dict_5p["C"])
                    umi_qual = Bc.get_seq_str(qual1, self.pattern_dict_5p["U"])
                    fh_5p_1.write(
                        "@{}\n{}\n+\n{}\n".format(header1, bc + umi, bc_qual + umi_qual)
                    )
            fh_5p_1.close()
            with pysam.FastxFile(fn2, persist=False) as fq:
                for entry2 in fq:
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    seq2 = utils.reverse_complement(seq2)
                    qual2 = qual2[::-1]
                    fh_5p_2.write("@{}\n{}\n+\n{}\n".format(header2, seq2, qual2))
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
        help="chemistry version",
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
