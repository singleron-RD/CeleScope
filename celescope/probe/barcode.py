import unittest
from collections import Counter, defaultdict
import sys

import pysam
import pandas as pd

import celescope.tools.parse_chemistry as parse_chemistry
from celescope.__init__ import HELP_DICT
from celescope.chemistry_dict import chemistry_dict
from celescope.tools import utils
from celescope.tools.step import Step, s_common


SEQUENCE_LENGTH = 150
PRIMER_OFFSET = 5


@utils.add_log
def fasta_to_dict(fasta, max_mismatch):
    seq_dict, _ = utils.read_fasta(fasta)
    mismatch_dict = {}
    length_dict = {}
    for seq_id, seq in seq_dict.items():
        mismatch_dict[seq_id] = parse_chemistry.create_mismatch_seqs(
            seq, max_mismatch=max_mismatch
        )
        length_dict[seq_id] = len(seq)

    return mismatch_dict, length_dict


@utils.add_log
def get_metrics_df(read_counter, umi_counter, total_reads, total_umis):
    all_probes = read_counter.keys()
    df = pd.DataFrame(index=sorted(all_probes))
    df["read"] = df.index.map(lambda x: read_counter.get(x, 0))
    df["umi"] = df.index.map(lambda x: umi_counter.get(x, 0))
    df["read_percent"] = df["read"] / total_reads * 100
    df["umi_percent"] = df["umi"] / total_umis * 100

    df["read_percent"] = df["read_percent"].round(2)
    df["umi_percent"] = df["umi_percent"].round(2)
    df = df.reset_index().rename(columns={"index": "id"})
    df = df.sort_values(by="read", ascending=False)
    return df


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
        self.bc_umi_length = sum(
            s.stop - s.start for s_list in self.pattern_dict.values() for s in s_list
        )
        sys.stderr.write(f"bc+umi+linker length: {self.bc_umi_length}\n")

        self.raw_list, self.mismatch_list = (
            parse_chemistry.create_mismatch_origin_dicts_from_whitelists(self.bc, 1)
        )
        # v3
        self.offset_runner = parse_chemistry.AutoRNA(self.fq1_list)
        #
        self.probe_mismatch_dict = {}
        if args.probe_fasta:
            self.probe_mismatch_dict, self.probe_length_dict = fasta_to_dict(
                args.probe_fasta, args.max_mismatch
            )

        self.primer_mismatch_dict = {}
        if args.primer_fasta:
            self.primer_mismatch_dict, self.primer_length_dict = fasta_to_dict(
                args.primer_fasta, args.max_mismatch
            )

        # output
        self.probe_count_file = f"{self.out_prefix}_probe_count.tsv"
        self.primer_count_file = f"{self.out_prefix}_primer_count.tsv"
        self.outs = [self.probe_count_file, self.primer_count_file]

    def get_bc_umi(self, seq):
        if self.chemistry == "GEXSCOPE-V3":
            offset = self.offset_runner.v3_offset(seq)
            seq = seq[offset:]
        elif self.chemistry == "flv_rna-V2":
            offset = self.offset_runner.flv_rna_v2_offset(seq)
            seq = seq[offset:]
        bc_list = [seq[x] for x in self.pattern_dict["C"]]
        if self.chemistry == "flv":
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
        all_umi_set = set()

        probe_read_counter = Counter()
        primer_read_counter = Counter()
        probe_umi_set = defaultdict(set)
        primer_umi_set = defaultdict(set)

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
                    all_umi_set.add((corrected_seq, umi))

                    for probe_id, seqs in self.probe_mismatch_dict.items():
                        seq_length = self.probe_length_dict[probe_id]
                        for start in range(
                            self.bc_umi_length, SEQUENCE_LENGTH - seq_length
                        ):
                            if e1.sequence[start : start + seq_length] in seqs:
                                probe_read_counter[probe_id] += 1
                                probe_umi_set[probe_id].add((corrected_seq, umi))
                                break

                    for primer_id, seqs in self.primer_mismatch_dict.items():
                        seq_length = self.primer_length_dict[primer_id]
                        for start in range(0, PRIMER_OFFSET):
                            if e2.sequence[start : start + seq_length] in seqs:
                                primer_read_counter[primer_id] += 1
                                primer_umi_set[primer_id].add((corrected_seq, umi))
                                break

                cb_quality_counter.update(
                    "".join([e1.quality[slice] for slice in self.pattern_dict["C"]])  # type: ignore
                )
                umi_quality_counter.update(
                    "".join([e1.quality[self.pattern_dict["U"][0]]])  # type: ignore
                )

        total_umi = len(all_umi_set)

        probe_umi_counter = {k: len(v) for k, v in probe_umi_set.items()}
        probe_df = get_metrics_df(
            probe_read_counter, probe_umi_counter, valid_reads, total_umi
        )
        probe_df.to_csv(self.probe_count_file, sep="\t", index=False)
        self.add_table(
            title="Probe",
            table_id="table_probe",
            df=probe_df,
        )

        primer_umi_counter = {k: len(v) for k, v in primer_umi_set.items()}
        primer_df = get_metrics_df(
            primer_read_counter, primer_umi_counter, valid_reads, total_umi
        )
        primer_df.to_csv(self.primer_count_file, sep="\t", index=False)
        self.add_table(
            title="Primer",
            table_id="table_primer",
            df=primer_df,
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
        self.add_metric(
            name="Total UMI",
            value=total_umi,
            help_info="total UMI count",
        )
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
    parser.add_argument("--probe_fasta", help="Probe fasta file.")
    parser.add_argument("--primer_fasta", help="Primer fasta file.")
    parser.add_argument(
        "--max_mismatch",
        help="Maximum mismatch in probe in primer sequence.",
        default=2,
        type=int,
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
