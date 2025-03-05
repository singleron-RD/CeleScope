"""Split fastq file according to barcode. Append UMI to read name with umi_separator semicolon ':'"""

import io
from celescope.tools.step import Step, s_common
from celescope.tools import utils
import celescope.tools.parse_chemistry as parse_chemistry
import pysam
from celescope.chemistry_dict import chemistry_dict
from celescope.__init__ import HELP_DICT
from celescope.bulk_rna.starsolo import get_barcode_sample


UMI_SEPARATOR = ":"


class Split_fastq(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")

        chemistry = parse_chemistry.get_chemistry(
            self.assay, args.chemistry, self.fq1_list
        )
        self.pattern_dict, bc = parse_chemistry.get_pattern_dict_and_bc(
            chemistry, args.pattern, args.whitelist
        )
        self.raw_list, self.mismatch_list = (
            parse_chemistry.create_mismatch_origin_dicts_from_whitelists(bc, 1)
        )

        self.barcode_sample = get_barcode_sample(bc[0], args.well_sample)

    @utils.add_log
    def split(self):
        self.fh_dict = {}
        for barcode, sample in self.barcode_sample.items():
            file_obj = utils.generic_open(
                f"{self.outdir}/{sample}.fq.gz", "wb", compresslevel=1
            )
            self.fh_dict[barcode] = io.BufferedWriter(
                file_obj, buffer_size=16 * 1024 * 1024
            )

        for fq1, fq2 in zip(self.fq1_list, self.fq2_list):
            fq1 = pysam.FastxFile(fq1)
            fq2 = pysam.FastxFile(fq2)
            for e1, e2 in zip(fq1, fq2):
                bc_list = [e1.sequence[x] for x in self.pattern_dict["C"]]
                valid, corrected, corrected_bc = parse_chemistry.check_seq_mismatch(
                    bc_list, self.raw_list, self.mismatch_list
                )
                if valid and corrected_bc in self.barcode_sample:
                    umi = e1.sequence[self.pattern_dict["U"][0]]
                    name = e1.name + UMI_SEPARATOR + umi
                    self.fh_dict[corrected_bc].write(
                        utils.fastq_line(name, e2.sequence, e2.quality).encode()
                    )

    def run(self):
        self.split()


def split_fastq(args):
    split_fastq_obj = Split_fastq(args)
    split_fastq_obj.run()


def get_opts_split_fastq(parser, sub_program):
    parser.add_argument(
        "--chemistry",
        help=HELP_DICT["chemistry"],
        choices=list(chemistry_dict.keys()),
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
        "--well_sample",
        help="tsv file of well numbers and sample names. The first column is well numbers and the second column is sample names.",
        required=True,
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
            "--barcodes",
            help="barcode file.",
            required=True,
        )
        s_common(parser)
    return parser
