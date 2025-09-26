import unittest
from collections import Counter

import pysam

import celescope.tools.parse_chemistry as parse_chemistry
from celescope.__init__ import HELP_DICT
from celescope.chemistry_dict import chemistry_dict
from celescope.tools import utils
from celescope.tools.step import Step, s_common

MAX_TRUST4_READ_PER_BARCODE = 80000


class Barcode(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        self.chemistry = parse_chemistry.get_chemistry(
            self.assay, args.chemistry, self.fq1_list
        )
        self.pattern_dict, self.bc = parse_chemistry.get_pattern_dict_and_bc(
            self.chemistry, args.pattern, args.whitelist
        )

        self.raw_list, self.mismatch_list = (
            parse_chemistry.create_mismatch_origin_dicts_from_whitelists(self.bc, 1)
        )
        # v3
        self.offset_runner = parse_chemistry.AutoRNA(self.fq1_list)
        # output
        self.out_fq2 = f"{self.outdir}/{self.sample}_2.fq"
        self.isflv = self.assay == "flv_trust4"
        if self.isflv:
            self.out_fq1 = f"{self.outdir}/{self.sample}_1.fq"
            self.match_barcodes = set(
                utils.get_barcode_from_match_dir(args.match_dir)[0]
            )

    def get_bc_umi(self, seq):
        if self.chemistry == "GEXSCOPE-V3":
            offset = self.offset_runner.v3_offset(seq)
            seq = seq[offset:]
        elif self.chemistry == "flv_rna-V2":
            offset = self.offset_runner.flv_rna_v2_offset(seq)
            seq = seq[offset:]
        bc_list = [seq[x] for x in self.pattern_dict["C"]]
        if self.isflv:
            bc_list = [utils.reverse_complement(bc) for bc in bc_list[::-1]]
        valid, corrected, corrected_seq = parse_chemistry.check_seq_mismatch(
            bc_list, self.raw_list, self.mismatch_list
        )
        if not valid:
            umi = None
        else:
            umi = seq[self.pattern_dict["U"][0]]
        return valid, corrected, corrected_seq, umi

    @utils.add_log
    def run(self):
        raw_reads = valid_reads = corrected_reads = 0
        # quality
        cb_quality_counter = Counter()
        umi_quality_counter = Counter()
        if self.isflv:
            out_fh1 = utils.generic_open(self.out_fq1, "wt")
            bc_read_counter = Counter()
            match_reads = 0

        with utils.generic_open(self.out_fq2, "wt") as out_fh2:
            for fq1, fq2 in zip(self.fq1_list, self.fq2_list):
                fq1 = pysam.FastxFile(fq1)
                fq2 = pysam.FastxFile(fq2)
                for e1, e2 in zip(fq1, fq2):
                    raw_reads += 1
                    valid, corrected, corrected_seq, umi = self.get_bc_umi(e1.sequence)
                    if valid:
                        valid_reads += 1
                        if corrected:
                            corrected_reads += 1
                        read_name = f"{corrected_seq}:{umi}:{raw_reads}"
                        if not self.isflv:
                            out_fh2.write(
                                utils.fastq_line(read_name, e2.sequence, e2.quality)
                            )  # type: ignore
                        else:
                            if corrected_seq not in self.match_barcodes:
                                continue
                            match_reads += 1
                            bc_read_counter[corrected_seq] += 1
                            seq1 = corrected_seq + umi
                            qual1 = "F" * len(seq1)
                            if (
                                bc_read_counter[corrected_seq]
                                <= MAX_TRUST4_READ_PER_BARCODE
                            ):
                                out_fh1.write(utils.fastq_line(read_name, seq1, qual1))
                                out_fh2.write(
                                    utils.fastq_line(read_name, e2.sequence, e2.quality)
                                )  # type: ignore

                    cb_quality_counter.update(
                        "".join([e1.quality[slice] for slice in self.pattern_dict["C"]])  # type: ignore
                    )
                    umi_quality_counter.update(
                        "".join([e1.quality[self.pattern_dict["U"][0]]])  # type: ignore
                    )

        self.add_metric(
            name="Raw Reads",
            value=raw_reads,
            help_info="total reads from FASTQ files",
        )
        self.add_metric(
            name="Valid Reads",
            value=valid_reads,
            total=raw_reads,
            help_info="reads with correct barcode",
        )
        self.add_metric(
            name="Corrected Reads",
            value=corrected_reads,
            total=raw_reads,
            help_info="Reads with barcodes that are not in the whitelist but are within one Hamming distance of it",
        )
        if not cb_quality_counter:
            q30_cb = q30_umi = 0
        else:
            q30_cb = sum(
                [
                    cb_quality_counter[k]
                    for k in cb_quality_counter
                    if Barcode.chr_to_int(k) >= 30
                ]
            ) / float(sum(cb_quality_counter.values()))
            q30_umi = sum(
                [
                    umi_quality_counter[k]
                    for k in umi_quality_counter
                    if Barcode.chr_to_int(k) >= 30
                ]
            ) / float(sum(umi_quality_counter.values()))
        self.add_metric(
            name="Q30 of Barcode",
            value=f"{round(q30_cb * 100,2)}%",
            help_info="percent of barcode base pairs with quality scores over Q30",
        )
        self.add_metric(
            name="Q30 of UMI",
            value=f"{round(q30_umi * 100, 2)}%",
            help_info="percent of UMI base pairs with quality scores over Q30",
        )

        if self.assay == "flv_trust4":
            self.add_metric(
                name="Valid Matched Reads",
                value=match_reads,
                total=valid_reads,
                help_info="reads match with flv_rna cell barcodes",
            )

            self.add_metric(
                name="Matched Barcodes",
                value=len(bc_read_counter),
                help_info="barcodes match with flv_rna",
            )

    @staticmethod
    def chr_to_int(chr, offset=33):
        """Convert Phred quality to int"""
        return ord(chr) - offset


@utils.add_log
def barcode(args):
    with Barcode(args, display_title="Demultiplexing") as runner:
        runner.run()


def get_opts_barcode(parser, sub_program=True):
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
        parser = s_common(parser)

    return parser


if __name__ == "__main__":
    unittest.main()
