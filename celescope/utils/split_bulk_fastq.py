"""Split bulk R1 and R2 fastq files by barcode for bulk_rna and bulk_vdj workflows."""

import io
import os
import sys

import pysam

import celescope.tools.parse_chemistry as parse_chemistry
from celescope.bulk_rna.starsolo import get_barcode_sample
from celescope.chemistry_dict import chemistry_dict
from celescope.__init__ import HELP_DICT
from celescope.tools import utils


class Split_bulk_fastq:
    """
    ## Features
    - Split R1 and R2 fastq files according to well barcodes.
    - Supports both `bulk_rna` and `bulk_vdj` workflows.

    ## Output
    - `{sample}_barcode_R1.fastq.gz` R1 fastq file for each well sample.
    - `{sample}_barcode_R2.fastq.gz` R2 fastq file for each well sample.
    """

    def __init__(self, args):
        self.args = args

    @utils.add_log
    def run(self):
        fq1_list = self.args.fq1.split(",")

        chemistry = parse_chemistry.get_chemistry(
            self.args.assay, self.args.chemistry, fq1_list
        )
        pattern_dict, bc = parse_chemistry.get_pattern_dict_and_bc(
            chemistry, self.args.pattern, self.args.whitelist
        )
        raw_list, mismatch_list = (
            parse_chemistry.create_mismatch_origin_dicts_from_whitelists(bc, 1)
        )
        barcode_sample = get_barcode_sample(bc[0], self.args.well_sample)

        os.makedirs(self.args.outdir, exist_ok=True)

        r1_fh_dict = {}
        r2_fh_dict = {}
        for barcode, sample in barcode_sample.items():
            r1_file = utils.generic_open(
                f"{self.args.outdir}/{sample}_{barcode}_R1.fastq.gz",
                "wb",
                compresslevel=1,
            )
            r2_file = utils.generic_open(
                f"{self.args.outdir}/{sample}_{barcode}_R2.fastq.gz",
                "wb",
                compresslevel=1,
            )
            r1_fh_dict[barcode] = io.BufferedWriter(
                r1_file, buffer_size=16 * 1024 * 1024
            )
            r2_fh_dict[barcode] = io.BufferedWriter(
                r2_file, buffer_size=16 * 1024 * 1024
            )

        valid_reads = 0
        total_reads = 0
        fq1_list = self.args.fq1.split(",")
        fq2_list = self.args.fq2.split(",")
        for fq1, fq2 in zip(fq1_list, fq2_list):
            fq1_obj = pysam.FastxFile(fq1)
            fq2_obj = pysam.FastxFile(fq2)
            for e1, e2 in zip(fq1_obj, fq2_obj):
                total_reads += 1
                if total_reads % 1000000 == 0:
                    sys.stderr.write(f"Processed {total_reads} reads\n")
                bc_list = [e1.sequence[x] for x in pattern_dict["C"]]
                valid, corrected, corrected_bc = parse_chemistry.check_seq_mismatch(
                    bc_list, raw_list, mismatch_list
                )
                if valid and corrected_bc in barcode_sample:
                    valid_reads += 1
                    r1_fh_dict[corrected_bc].write(
                        utils.fastq_line(e1.name, e1.sequence, e1.quality).encode()
                    )
                    r2_fh_dict[corrected_bc].write(
                        utils.fastq_line(e2.name, e2.sequence, e2.quality).encode()
                    )

        for fh in r1_fh_dict.values():
            fh.close()
        for fh in r2_fh_dict.values():
            fh.close()

        sys.stderr.write(
            f"Split {valid_reads}/{total_reads} valid reads to {len(barcode_sample)} samples\n"
        )


@utils.add_log
def split_bulk_fastq(args):
    Split_bulk_fastq(args).run()


def get_opts_split_bulk_fastq(parser, sub_program=True):
    parser.add_argument(
        "--assay",
        help="Assay type. Determines how to parse the barcode from R1 reads.",
        choices=["bulk_rna", "bulk_vdj"],
        required=True,
    )
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
            "--outdir",
            help="Output directory.",
            required=True,
        )
    return parser
