"""
bulk-seq vdj mapping
"""

import pandas as pd
import subprocess

from celescope.tools import utils
from celescope.vdj import mapping_vdj as super_vdj


class Mapping_vdj(super_vdj.Mapping_vdj):
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

        # out
        self.airr_out = f"{self.out_prefix}_airr.tsv"
        self.productive_file = f"{self.out_prefix}_productive.tsv"

    @utils.add_log
    def igblast(self, chain, ig_seqtype):
        cmd = (
            f"igblastn -query {self.args.fasta} "
            f"-organism {self.species} "
            f"-ig_seqtype {ig_seqtype} "
            f"-auxiliary_data optional_file/{self.species}_gl.aux "
            f"-num_threads {self.args.thread} "
            f"-germline_db_V {self.ref_path}/{chain}V.fa "
            f"-germline_db_D {self.ref_path}/{chain}D.fa "
            f"-germline_db_J {self.ref_path}/{chain}J.fa "
            "-domain_system imgt -show_translation -outfmt 19 "  # outfmt19 is an AIRR tab-delimited file, IgBLAST v1.9.0 or higher required.
            f"-out {self.airr_out} "
        )

        self.igblast.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def mapping_summary(self):
        df = pd.read_csv(self.airr_out, sep="\t")
        df.fillna("", inplace=True)
        total_reads = df.shape[0]

        self.add_metric(name="Species", value=self.species, help_info="Human or Mouse")

        df = df[(df["v_call"] != "") | (df["d_call"] != "") | (df["j_call"] != "")]
        self.add_metric(
            name="UMIs Mapped to Any VDJ Gene",
            value=df.shape[0],
            total=total_reads,
            help_info="UMIs Mapped to any germline VDJ gene segments",
        )

        # Reads with CDR3
        df_cdr3 = df[(df["cdr3_aa"] != "") & (df["junction_aa"] != "")]
        self.add_metric(
            name="UMIs with CDR3",
            value=df_cdr3.shape[0],
            total=total_reads,
            help_info="UMIs with CDR3 sequence",
        )

        # Reads with Correct CDR3
        df_correct_cdr3 = df_cdr3[
            ~(df_cdr3["cdr3_aa"].str.contains(r"\*"))
            & ~(df_cdr3["cdr3_aa"].str.contains("X"))
        ]
        self.add_metric(
            name="UMIs with Correct CDR3",
            value=df_correct_cdr3.shape[0],
            total=total_reads,
            help_info="UMIs with CDR3 might have stop codon and these Reads are classified as incorrect",
        )

        # Reads Mapped Confidently To VJ Gene
        df_confident = df_correct_cdr3[df_correct_cdr3["productive"] == "T"]
        self.add_metric(
            name="UMIs Mapped Confidently to VJ Gene",
            value=df_confident.shape[0],
            total=total_reads,
            help_info="UMIs with productive rearrangement mapped to VJ gene pairs and with correct CDR3",
        )

        # Reads Mapped Confidently to each chain
        for chain in self.chains:
            df_chain = df_confident[df_confident.locus == chain]
            self.add_metric(
                name=f"UMIs Mapped to {chain}",
                value=df_chain.shape[0],
                total=total_reads,
                help_info=f"UMIs mapped confidently to {chain}",
            )

        # output file
        df_confident["barcode"] = df_confident["sequence_id"].apply(
            lambda x: x.split(":")[0]
        )
        df_VJ = df_confident[
            [
                "barcode",
                "sequence_id",
                "locus",
                "v_call",
                "d_call",
                "j_call",
                "junction",
                "junction_aa",
            ]
        ]
        df_VJ.rename(
            columns={
                "locus": "chain",
                "v_call": "bestVGene",
                "d_call": "bestDGene",
                "j_call": "bestJGene",
                "junction": "nSeqCDR3",
                "junction_aa": "aaSeqCDR3",
            },
            inplace=True,
        )

        for i in ["bestVGene", "bestDGene", "bestJGene"]:
            df_VJ[i] = df_VJ[i].apply(lambda x: x.split("*")[0])

        # CDR3 sequence have at least 5 amino acids and start with C
        df_VJ = df_VJ[df_VJ["aaSeqCDR3"].str.len() > 5]
        df_VJ = df_VJ[df_VJ["aaSeqCDR3"].str.startswith("C")]
        df_VJ.to_csv(self.productive_file, sep="\t", index=False)

    def run(self):
        # run igblstn
        if self.seqtype == "TCR":
            self.igblast(chain="TR", ig_seqtype="TCR")
        else:
            self.igblast(chain="IG", ig_seqtype="Ig")
        self.mapping_summary()


def mapping_vdj(args):
    with Mapping_vdj(args, display_title="Mapping") as runner:
        runner.run()


def get_opts_mapping_vdj(parser, sub_program):
    super_vdj.get_opts_mapping_vdj(parser, sub_program)
