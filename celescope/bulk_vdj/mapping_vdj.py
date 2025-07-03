"""
bulk-seq vdj mapping
"""

import subprocess
import csv
import os

from celescope.tools import utils, step
from celescope.vdj import mapping_vdj as super_vdj
from celescope.vdj.__init__ import CHAINS
from celescope.bulk_rna.starsolo import get_barcode_sample, get_well_barcode
import celescope.tools.parse_chemistry as parse_chemistry
from celescope.__init__ import HELP_DICT
from celescope.chemistry_dict import chemistry_dict
import pandas as pd


def simpson_diversity(data):
    """Given a hash { 'species': count } , returns the Simpson Diversity Index
    >>> simpson_di({'a': 10, 'b': 20, 'c': 30,})
    0.3888888888888889
    """

    def p(n, N):
        """Relative abundance"""
        if n == 0:
            return 0
        else:
            return float(n) / N

    N = sum(data.values())

    return sum(p(n, N) ** 2 for n in data.values() if n != 0)


def inverse_simpson_diversity(data):
    """Given a hash { 'species': count } , returns the inverse Simpson Diversity Index
    >>> inverse_simpson_diversity({'a': 10, 'b': 20, 'c': 30,})
    2.571428571428571
    """
    return float(1) / simpson_diversity(data)


class Mapping_vdj(step.Step):
    """
    ## Features
    - Align R2 reads to IMGT(http://www.imgt.org/) database sequences with blast.
    ## Output
    - `{sample}_airr.tsv` The alignment result of each read.
    A tab-delimited file compliant with the AIRR Rearrangement schema(https://docs.airr-community.org/en/stable/datarep/rearrangements.html)
    - `{sample}_produtive.tsv` Including all productive chains.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.chains = CHAINS[args.type]

        chemistry = self.get_slot_key(
            slot="metrics", step_name="sample", key="Chemistry"
        )

        self.pattern_dict, bc = parse_chemistry.get_pattern_dict_and_bc(
            chemistry, args.pattern, args.whitelist
        )
        self.barcode_sample = get_barcode_sample(bc[0], args.well_sample)
        well_barcode = get_well_barcode(bc[0])
        self.barcode_well = {v: k for k, v in well_barcode.items()}

        # out
        self.productive_dir = f"{self.outdir}/productive"
        if not os.path.exists(self.productive_dir):
            os.makedirs(self.productive_dir)
        self.annotation_dir = f"{self.outdir}/annotation"
        if not os.path.exists(self.annotation_dir):
            os.makedirs(self.annotation_dir)
        self.clonotypes_dir = f"{self.outdir}/clonotypes"
        if not os.path.exists(self.clonotypes_dir):
            os.makedirs(self.clonotypes_dir)
        self.airr_out = f"{self.out_prefix}_airr.tsv"
        self.annotation_csv = f"{self.out_prefix}_filtered_annotations.csv"
        self.clonotypes_csv = f"{self.out_prefix}_clonotypes.csv"
        self.outs = [self.annotation_csv, self.clonotypes_csv]

    @utils.add_log
    def igblast(self):
        if self.args.type == "TCR":
            chain = "TR"
            ig_seqtype = "TCR"
        elif self.args.type == "BCR":
            chain = "IG"
            ig_seqtype = "Ig"
        cmd = (
            f"igblastn -query {self.args.fasta} "
            f"-organism {self.args.species} "
            f"-ig_seqtype {ig_seqtype} "
            f"-auxiliary_data optional_file/{self.args.species}_gl.aux "
            f"-num_threads {self.args.thread} "
            f"-germline_db_V {self.args.ref_path}/{chain}V.fa "
            f"-germline_db_D {self.args.ref_path}/{chain}D.fa "
            f"-germline_db_J {self.args.ref_path}/{chain}J.fa "
            "-domain_system imgt -show_translation -outfmt 19 "  # outfmt19 is an AIRR tab-delimited file, IgBLAST v1.9.0 or higher required.
            f"-out {self.airr_out} "
        )
        self.igblast.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def process_airr(self):
        """
        read airr line by line, collect metrics, write productive to seperate wells
        """

        tsv_handles = {
            barcode: open(f"{self.productive_dir}/{sample}_productive.tsv", "wt")
            for barcode, sample in self.barcode_sample.items()
        }
        tsv_headers = [
            "barcode",
            "well",
            "sample",
            "umi",
            "chain",
            "bestVGene",
            "bestDGene",
            "bestJGene",
            "nSeqCDR3",
            "aaSeqCDR3",
        ]
        for f in tsv_handles.values():
            f.write("\t".join(tsv_headers) + "\n")

        consensus_metrics_df = pd.read_csv(
            self.args.consensus_metrics_file, sep="\t", index_col=0
        )
        consensus_metrics_dict = consensus_metrics_df.to_dict(orient="index")
        metrics = utils.nested_defaultdict(dim=2)
        for barcode in consensus_metrics_dict:
            if barcode in self.barcode_sample:
                metrics[barcode]["well"] = self.barcode_well[barcode]
                metrics[barcode]["sample"] = self.barcode_sample[barcode]
            else:
                metrics[barcode]["well"] = ""
                metrics[barcode]["sample"] = ""
            metrics[barcode]["n_clonotypes"] = 0
            metrics[barcode]["diversity"] = 0
            metrics[barcode]["read"] = consensus_metrics_dict[barcode]["n_read"]
            metrics[barcode]["umi"] = 0
            metrics[barcode]["umi_mapped"] = 0
            metrics[barcode]["umi_confident"] = 0
            for chain in self.chains:
                metrics[barcode][f"umi_confident_{chain}"] = 0

        total_umi = 0
        with open(self.airr_out, "rt") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                total_umi += 1
                barcode, umi = row["sequence_id"].split(":")[0:2]
                metrics[barcode]["umi"] += 1
                if row["v_call"] != "" or row["d_call"] != "" or row["j_call"] != "":
                    metrics[barcode]["umi_mapped"] += 1
                else:
                    continue
                if (
                    row["productive"] == "T"
                    and row["junction"] != ""
                    and "N" not in row["junction"]
                ):
                    metrics[barcode]["umi_confident"] += 1
                    if row["locus"] in self.chains:
                        metrics[barcode][f"umi_confident_{row['locus']}"] += 1
                else:
                    continue
                if barcode in self.barcode_sample:
                    for gene in ["v_call", "d_call", "j_call"]:
                        row[gene] = row[gene].split("*")[0]
                    line = [
                        barcode,
                        str(self.barcode_well[barcode]),
                        self.barcode_sample[barcode],
                        umi,
                        row["locus"],
                        row["v_call"],
                        row["d_call"],
                        row["j_call"],
                        row["junction"],
                        row["junction_aa"],
                    ]
                    tsv_handles[barcode].write("\t".join(line) + "\n")

        umi_mapped = sum(v for barcode, k in metrics.items() for v in [k["umi_mapped"]])
        self.add_metric(
            name="UMIs Mapped to Any VDJ Gene",
            value=umi_mapped,
            total=total_umi,
            help_info="UMIs Mapped to any germline VDJ gene segments",
        )

        umi_confident = sum(
            v for barcode, k in metrics.items() for v in [k["umi_confident"]]
        )
        self.add_metric(
            name="UMIs Mapped Confidently to VJ Gene",
            value=umi_confident,
            total=total_umi,
            help_info="UMIs with productive rearrangement mapped to VJ gene pairs and without N in the junction sequence",
        )

        for chain in self.chains:
            umi_confident_chain = sum(
                v
                for barcode, k in metrics.items()
                for v in [k[f"umi_confident_{chain}"]]
            )
            self.add_metric(
                name=f"UMIs Mapped Confidently to {chain}",
                value=umi_confident_chain,
                total=total_umi,
            )

        self.metrics = metrics

    @staticmethod
    def create_well_annotation(productive_file, out_file):
        df = pd.read_csv(productive_file, sep="\t")
        df.fillna("None", inplace=True)
        df = (
            df.groupby(
                [
                    "barcode",
                    "well",
                    "sample",
                    "chain",
                    "bestVGene",
                    "bestDGene",
                    "bestJGene",
                    "nSeqCDR3",
                    "aaSeqCDR3",
                ]
            )
            .size()
            .reset_index(name="umis")
        )
        out_df = df.assign(
            is_cell=True,
            high_confidence=True,
            c_gene="None",
            full_length=True,
            productive=True,
        ).rename(
            columns={
                "bestVGene": "v_gene",
                "bestDGene": "d_gene",
                "bestJGene": "j_gene",
                "aaSeqCDR3": "cdr3",
                "nSeqCDR3": "cdr3_nt",
            }
        )[
            [
                "barcode",
                "well",
                "sample",
                "is_cell",
                "high_confidence",
                "chain",
                "v_gene",
                "d_gene",
                "j_gene",
                "c_gene",
                "full_length",
                "productive",
                "cdr3",
                "cdr3_nt",
                "umis",
            ]
        ]

        out_df.sort_values(by=["umis"], ascending=False, inplace=True)
        out_df.to_csv(out_file, sep=",", index=False)

    @utils.add_log
    def create_annotation(self):
        for barcode, sample in self.barcode_sample.items():
            productive_file = f"{self.productive_dir}/{sample}_productive.tsv"
            out_file = f"{self.annotation_dir}/{sample}_annotation.csv"
            self.create_well_annotation(productive_file, out_file)

        annotation_files = [
            f"{self.annotation_dir}/{sample}_annotation.csv"
            for sample in self.barcode_sample.values()
        ]
        utils.merge_table_files(annotation_files, self.annotation_csv)

    @staticmethod
    def create_well_clonotypes(annotation_file, out_file):
        df = pd.read_csv(annotation_file)
        df = (
            df.groupby(["barcode", "well", "sample", "chain", "cdr3"])["umis"]
            .sum()
            .reset_index(name="umis")
        )
        total_umis = sum(df["umis"])
        df["percent"] = df["umis"].apply(lambda x: f"{x / total_umis * 100:.2f}%")
        df.sort_values(by=["umis"], ascending=False, inplace=True)
        df.to_csv(out_file, index=False)
        diversity = inverse_simpson_diversity(dict(zip(df["cdr3"], df["umis"])))
        n_clonotypes = len(df)
        return n_clonotypes, diversity

    @utils.add_log
    def create_clonotypes(self):
        for barcode, sample in self.barcode_sample.items():
            annotation_file = f"{self.annotation_dir}/{sample}_annotation.csv"
            out_file = f"{self.clonotypes_dir}/{sample}_clonotypes.csv"
            n_clonotypes, diversity = self.create_well_clonotypes(
                annotation_file, out_file
            )
            self.metrics[barcode]["n_clonotypes"] = n_clonotypes
            self.metrics[barcode]["diversity"] = round(diversity, 2)

        clonotypes_files = [
            f"{self.clonotypes_dir}/{sample}_clonotypes.csv"
            for sample in self.barcode_sample.values()
        ]
        utils.merge_table_files(clonotypes_files, self.clonotypes_csv)

        df = pd.read_csv(self.clonotypes_csv, sep=",")
        df = df.drop(["barcode"], axis=1)
        df = df.groupby("sample", group_keys=False).apply(
            lambda x: x.sort_values("umis", ascending=False).head(10)
        )
        self.add_table(
            title="Top 10 Clonotypes",
            table_id="clonotypes",
            df=df,
        )

    @utils.add_log
    def add_metrics(self):
        metrics_df = pd.DataFrame.from_dict(self.metrics, orient="index")
        # add percent
        umi_cols = [
            "umi_mapped",
            "umi_confident",
        ]
        umi_cols.extend([f"umi_confident_{chain}" for chain in self.chains])
        for col in umi_cols:
            metrics_df[col] = metrics_df.apply(
                lambda row: f"{row[col]}({round(row[col]/row['umi']*100, 2)}%)", axis=1
            )

        metrics_df.sort_values("read", inplace=True)
        metrics_df.to_csv(f"{self.out_prefix}_well_metrics.tsv", sep="\t")
        well_metrics_df = metrics_df[
            metrics_df["sample"].isin(self.barcode_sample.values())
        ]
        well_reads = well_metrics_df["read"].sum()
        total_reads = metrics_df["read"].sum()

        self.add_metric(
            name="Fraction of Reads in Wells",
            value=well_reads / total_reads,
            value_type="fraction",
            help_info="Lower value indicates wrong well_sample file or high ambient contamination.",
        )
        well_metrics_df.sort_values("well", inplace=True)

        cols = ["well", "sample", "n_clonotypes", "diversity"]
        well_metrics_df = well_metrics_df[
            cols + [col for col in well_metrics_df.columns if col not in cols]
        ]
        self.add_table(
            title="Well Metrics",
            table_id="well_metrics",
            df=well_metrics_df,
            help="Diversity: inverse Simpson diversity index",
        )

    def run(self):
        self.igblast()
        self.process_airr()
        self.create_annotation()
        self.create_clonotypes()
        self.add_metrics()


def mapping_vdj(args):
    with Mapping_vdj(args, display_title="Mapping") as runner:
        runner.run()


def get_opts_mapping_vdj(parser, sub_program):
    parser.add_argument(
        "--well_sample",
        help="tsv file of well numbers and sample names. The first column is well numbers and the second column is sample names.",
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
    super_vdj.get_opts_mapping_vdj(parser, sub_program)
    if sub_program:
        parser.add_argument("--consensus_metrics_file", required=True)
