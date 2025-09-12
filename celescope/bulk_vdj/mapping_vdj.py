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
import numpy as np


def cal_chao1(counts):
    F1 = sum(1 for c in counts if c == 1)
    F2 = sum(1 for c in counts if c == 2)
    S_obs = len([c for c in counts if c > 0])

    if F2 == 0:
        return S_obs + F1 * (F1 - 1) / (2 * (F2 + 1))
    else:
        return S_obs + F1**2 / (2 * F2)


def cal_inverse_simpson(counts):
    """
    >>> inverse_simpson([10,20,30])
    2.57
    """
    counts = np.array(counts)
    if np.sum(counts) == 0:
        return 0.0
    proportions = counts / np.sum(counts)
    return round(1.0 / np.sum(proportions**2), 2)


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
        self.raw_dir = f"{self.outdir}/raw"
        if not os.path.exists(self.raw_dir):
            os.makedirs(self.raw_dir)
        self.annotation_dir = f"{self.outdir}/annotation"
        if not os.path.exists(self.annotation_dir):
            os.makedirs(self.annotation_dir)
        self.clonotypes_dir = f"{self.outdir}/clonotypes"
        if not os.path.exists(self.clonotypes_dir):
            os.makedirs(self.clonotypes_dir)
        self.airr_out = f"{self.out_prefix}_airr.tsv"
        self.annotation_csv = f"{self.out_prefix}_filtered_annotations.csv"
        self.clonotypes_csv = f"{self.out_prefix}_clonotypes.csv"
        self.outs = [
            self.annotation_csv,
            self.clonotypes_csv,
            self.annotation_dir,
            self.clonotypes_dir,
        ]

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

        raw_handles = {
            barcode: open(f"{self.raw_dir}/{sample}_raw.csv", "wt")
            for barcode, sample in self.barcode_sample.items()
        }
        tsv_headers = [
            "barcode",
            "well",
            "sample",
            "umi",
            "chain",
            "v_gene",
            "d_gene",
            "j_gene",
            "c_gene",
            "productive",
            "cdr3",
            "cdr3_nt",
        ]
        for f in raw_handles.values():
            f.write(",".join(tsv_headers) + "\n")

        consensus_metrics_df = pd.read_csv(
            self.args.consensus_metrics_file, sep="\t", index_col=0
        )
        consensus_metrics_dict = consensus_metrics_df.to_dict(orient="index")
        metrics = utils.nested_defaultdict(dim=2)
        for barcode in consensus_metrics_dict:
            metrics[barcode]["well"] = self.barcode_well[barcode]
            if barcode in self.barcode_sample:
                metrics[barcode]["sample"] = self.barcode_sample[barcode]
            else:
                metrics[barcode]["sample"] = ""
            metrics[barcode]["n_clonotypes"] = 0
            metrics[barcode]["inverse_simpson"] = 0
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
                    and len(row["junction_aa"]) > 5
                    and row["junction_aa"][0] == "C"
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
                        "None",
                        "True",
                        row["junction_aa"],
                        row["junction"],
                    ]
                    raw_handles[barcode].write(",".join(line) + "\n")

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
            help_info="UMIs with productive rearrangement mapped to VJ gene pairs and with correct juntion sequence(no N in nt, aa length > 5, aa starts with C)",
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
    def create_well_annotation_clono(raw_file, annotation_file, clono_file):
        df = pd.read_csv(raw_file)
        df = (
            df.groupby(
                [
                    "barcode",
                    "well",
                    "sample",
                    "chain",
                    "v_gene",
                    "d_gene",
                    "j_gene",
                    "productive",
                    "cdr3_nt",
                    "cdr3",
                ],
                dropna=False,
            )
            .size()
            .reset_index(name="umis")
        )
        sample = df["sample"][0]
        df.sort_values("umis", ascending=False, inplace=True)
        df["raw_clonotype_id"] = ["{}_{}".format(sample, i + 1) for i in range(len(df))]
        df["barcode"] = [
            "{}_{}".format(barcode, i + 1)
            for barcode, i in zip(df["barcode"], range(len(df)))
        ]
        total_umis = sum(df["umis"])
        df["percent"] = df["umis"].apply(lambda x: f"{x / total_umis * 100:.2f}%")
        df.to_csv(annotation_file, index=False)

        df_clono = df[
            [
                "well",
                "sample",
                "chain",
                "v_gene",
                "d_gene",
                "j_gene",
                "cdr3_nt",
                "cdr3",
                "raw_clonotype_id",
                "umis",
                "percent",
            ]
        ]
        df_clono.to_csv(clono_file, index=False)
        inverse_simpson = cal_inverse_simpson(df_clono["umis"])
        n_clonotypes = len(df_clono)
        return n_clonotypes, inverse_simpson

    @utils.add_log
    def create_clonotypes(self):
        for barcode, sample in self.barcode_sample.items():
            raw_file = f"{self.raw_dir}/{sample}_raw.csv"
            annotation_file = f"{self.annotation_dir}/{sample}_annotation.csv"
            clono_file = f"{self.clonotypes_dir}/{sample}_clonotypes.csv"
            n_clonotypes, inverse_simpson = self.create_well_annotation_clono(
                raw_file, annotation_file, clono_file
            )
            self.metrics[barcode]["n_clonotypes"] = n_clonotypes
            self.metrics[barcode]["inverse_simpson"] = round(inverse_simpson, 2)

    @utils.add_log
    def merge_files(self):
        annotation_files = [
            f"{self.annotation_dir}/{sample}_annotation.csv"
            for sample in self.barcode_sample.values()
        ]
        utils.merge_table_files(annotation_files, self.annotation_csv)

        clonotypes_files = [
            f"{self.clonotypes_dir}/{sample}_clonotypes.csv"
            for sample in self.barcode_sample.values()
        ]
        utils.merge_table_files(clonotypes_files, self.clonotypes_csv)

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
                lambda row: (
                    f"{row[col]}(0.0%)"
                    if row["umi"] == 0
                    else f"{row[col]}({round(row[col] / row['umi'] * 100, 2)}%)"
                ),
                axis=1,
            )

        metrics_df.sort_values("read", ascending=False, inplace=True)
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

        cols = ["well", "sample", "n_clonotypes", "inverse_simpson"]
        well_metrics_df = well_metrics_df[
            cols + [col for col in well_metrics_df.columns if col not in cols]
        ]
        self.add_table(
            title="Well Metrics",
            table_id="well_metrics",
            df=well_metrics_df,
            help="n_clonotypes: total number of clonotypes<br>inverse_simpson: inverse Simpson diversity index<br>",
        )

        df = pd.read_csv(self.clonotypes_csv, sep=",")
        df = df.groupby("sample", group_keys=False).head(10)
        self.add_table(
            title="Top 10 Clonotypes",
            table_id="clonotypes",
            df=df,
        )

    def run(self):
        self.igblast()
        self.process_airr()
        self.create_clonotypes()
        self.merge_files()
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
