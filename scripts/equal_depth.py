"""keep the same amount of reads in each well barcode"""

import argparse
from collections import defaultdict

import pysam

import celescope.tools.parse_chemistry as parse_chemistry
from celescope.__init__ import HELP_DICT
from celescope.chemistry_dict import chemistry_dict
from celescope.tools import utils


def main():
    parser = get_parser(argparse.ArgumentParser())
    args = parser.parse_args()

    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")

    chemistry = parse_chemistry.get_chemistry("bulk_rna", args.chemistry, fq1_list)
    pattern_dict, bc = parse_chemistry.get_pattern_dict_and_bc(
        chemistry, args.pattern, args.whitelist
    )
    raw_list, mismatch_list = (
        parse_chemistry.create_mismatch_origin_dicts_from_whitelists(bc, 1)
    )

    barcode_readcount = defaultdict(int)
    out_r1_fn = f"{args.outdir}/out_R1.fq.gz"
    out_r2_fn = f"{args.outdir}/out_R2.fq.gz"
    out_r1 = utils.generic_open(out_r1_fn, "wt", compresslevel=1)
    out_r2 = utils.generic_open(out_r2_fn, "wt", compresslevel=1)

    for fq1, fq2 in zip(fq1_list, fq2_list):
        fq1 = pysam.FastxFile(fq1)
        fq2 = pysam.FastxFile(fq2)
        for e1, e2 in zip(fq1, fq2):
            bc_list = [e1.sequence[x] for x in pattern_dict["C"]]
            valid, corrected, corrected_bc = parse_chemistry.check_seq_mismatch(
                bc_list, raw_list, mismatch_list
            )
            if valid:
                barcode_readcount[corrected_bc] += 1
                if barcode_readcount[corrected_bc] <= args.max_read:
                    out_r1.write(utils.fastq_line(e1.name, e1.sequence, e1.quality))
                    out_r2.write(utils.fastq_line(e2.name, e2.sequence, e2.quality))

    out_r1.close()
    out_r2.close()


def get_parser(parser):
    parser = argparse.ArgumentParser()
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
        "--max_read",
        help="keep maximum this number of reads per well in the output fastq files",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--outdir",
        help="Output directory",
        default="./",
    )
    return parser


if __name__ == "__main__":
    main()
