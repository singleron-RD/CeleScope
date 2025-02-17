import os
import glob
import sys
import unittest


from celescope.tools import utils
from celescope.tools.__init__ import PATTERN_DICT
from celescope.__init__ import ROOT_PATH, HELP_DICT
from celescope.tools.step import Step, s_common

from sccore.parse_protocol import AutoRNA, parse_pattern
from sccore.extract import Extract


def get_whitelist(chemistry):
    """
    returns: [bclists]
    """
    pattern = PATTERN_DICT[chemistry]
    pattern_dict = Barcode.parse_pattern(pattern)
    repeat = len(pattern_dict["C"])
    root_dir = f"{ROOT_PATH}/data/chemistry/{chemistry}"
    bclist = f"{root_dir}/bclist"
    bclist1 = f"{root_dir}/bclist1"
    if os.path.exists(bclist1):
        return [f"{root_dir}/bclist{i}" for i in range(1, repeat + 1)]
    elif os.path.exists(bclist):
        return [bclist] * repeat
    else:
        sys.exit(f"No bclist found for chemistry {chemistry} under {root_dir}")


def get_linker_whitelist_file(chemistry, wells=""):
    """Return (linker file path, whitelist file path)"""
    try:
        if chemistry == "bulk_rna":
            linker_f = None
            whitelist_f = f"{ROOT_PATH}/data/chemistry/{chemistry}/bclist{wells}"
        else:
            linker_f = glob.glob(f"{ROOT_PATH}/data/chemistry/{chemistry}/linker*")[0]
            whitelist_f = f"{ROOT_PATH}/data/chemistry/{chemistry}/bclist"
    except IndexError:
        return None, None
    return linker_f, whitelist_f


def get_protocol_and_dict(assay, fq1_list, chemistry):
    if assay in ["bulk_rna", "bulk_vdj"]:
        return assay


class Barcode(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        if args.chemistry == "auto":
            self.chemistry, self.protocol_dict = AutoRNA(self.fq1_list).run()
        elif args.chemistry == "customized":
            self.chemistry = "customized"
            self.protocol_dict = {
                "pattern_dict": parse_pattern(args.pattern),
                "bc": args.whitelist.split(","),
            }

    def run(self):
        runner = Extract(
            self.fq1_list,
            self.fq2_list,
            self.sample,
            self.protocol_dict,
            self.chemistry,
            self.outdir,
        )
        runner.run()


@utils.add_log
def barcode(args):
    with Barcode(args, display_title="Demultiplexing") as runner:
        runner.run()


def get_opts_barcode(parser, sub_program=True):
    parser.add_argument(
        "--chemistry",
        help="Predefined (pattern, barcode whitelist, linker whitelist) combinations. "
        + HELP_DICT["chemistry"],
        choices=list(PATTERN_DICT.keys()),
        default="auto",
    )
    parser.add_argument(
        "--pattern",
        help="""The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences)  
- `U`: UMI    
- `T`: poly T""",
    )
    parser.add_argument(
        "--whitelist",
        help="Cell barcode whitelist file path, one cell barcode per line.",
    )
    parser.add_argument(
        "--wells", help="The AccuraCode wells used (384 or 96).", type=int, default=384
    )
    if sub_program:
        parser.add_argument(
            "--fq1",
            help="R1 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser.add_argument(
            "--fq2",
            help="R2 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser.add_argument(
            "--match_dir", help="Matched scRNA-seq directory, required for flv_trust4"
        )
        parser.add_argument(
            "--stdout", help="Write output to standard output", action="store_true"
        )
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
