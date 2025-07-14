"""Split fastq file according to barcode. Append UMI to read name with umi_separator semicolon ':'"""

import io
from celescope.tools.step import Step, s_common
from celescope.tools import utils
import celescope.tools.parse_chemistry as parse_chemistry
import pysam
import sys
from celescope.chemistry_dict import chemistry_dict
from celescope.__init__ import HELP_DICT
from celescope.bulk_rna.starsolo import get_barcode_sample


UMI_SEPARATOR = ":"


class Split(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.split_fastq = args.split_fastq
        self.split_bam = args.split_bam

        fq1_list = args.fq1.split(",")
        chemistry = parse_chemistry.get_chemistry(self.assay, args.chemistry, fq1_list)
        self.pattern_dict, bc = parse_chemistry.get_pattern_dict_and_bc(
            chemistry, args.pattern, args.whitelist
        )
        self.raw_list, self.mismatch_list = (
            parse_chemistry.create_mismatch_origin_dicts_from_whitelists(bc, 1)
        )

        self.barcode_sample = get_barcode_sample(bc[0], args.well_sample)

    @utils.add_log
    def run(self):
        fq_dict = {}
        bam_dict = {}
        all_bam = pysam.AlignmentFile(self.args.bam, "rb")
        for barcode, sample in self.barcode_sample.items():
            if self.split_fastq:
                file_obj = utils.generic_open(
                    f"{self.outdir}/{sample}.fq.gz", "wb", compresslevel=1
                )
                fq_dict[barcode] = io.BufferedWriter(
                    file_obj, buffer_size=16 * 1024 * 1024
                )
            if self.split_bam:
                bam_dict[barcode] = pysam.AlignmentFile(
                    f"{self.outdir}/{sample}.bam", "wb", header=all_bam.header
                )

        dup_align_read_names = set()
        for segment in all_bam:
            cb = segment.get_tag("CB")
            if cb == "-":
                continue
            if cb not in self.barcode_sample:
                continue
            if self.split_bam:
                bam_dict[cb].write(segment)
            if segment.get_tag("NH") > 1:
                if segment.query_name in dup_align_read_names:
                    continue
                else:
                    dup_align_read_names.add(segment.query_name)
            if self.split_fastq:
                umi = segment.get_tag("UB")
                name = segment.query_name + UMI_SEPARATOR + umi
                seq = segment.get_forward_sequence()
                qual = "".join(chr(q + 33) for q in segment.get_forward_qualities())
                fq_dict[cb].write(utils.fastq_line(name, seq, qual).encode())


@utils.add_log
def split(args):
    if not args.split_fastq:
        sys.stderr.write("--split_fastq not set. Skip split_fastq\n")
    if not args.split_bam:
        sys.stderr.write("--split_bam not set. Skip split_bam\n")
    if not (args.split_fastq or args.split_bam):
        return
    Split(args).run()


def get_opts_split(parser, sub_program):
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
    parser.add_argument(
        "--split_fastq",
        help="Split fastq file according to well barcodes. Append UMI to read name with umi_separator semicolon ':'",
        action="store_true",
    )
    parser.add_argument(
        "--split_bam",
        help="Split BAM file according to well barcodes.",
        action="store_true",
    )
    if sub_program:
        parser.add_argument(
            "--fq1",
            help="R1 fastq file. Multiple files are separated by comma.",
            required=True,
        )
        parser.add_argument(
            "--bam",
            help="BAM file",
            required=True,
        )
        parser.add_argument(
            "--barcodes",
            help="barcode file.",
            required=True,
        )
        s_common(parser)
    return parser
