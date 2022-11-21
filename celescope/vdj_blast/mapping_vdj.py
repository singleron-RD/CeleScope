'''
vdj mapping
'''
import pandas as pd
import pysam
import os
import subprocess

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj_blast.__init__ import CHAINS


class Mapping_vdj(Step):

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.soft_path = args.soft_path
        self.ref_path = args.ref_path
        self.seqtype = args.type
        self.species = args.species
        self.chains = CHAINS[self.seqtype]
        self.fasta = os.path.abspath(args.fasta)

        # out
        self.airr_out = os.path.abspath(f"{self.out_prefix}_airr.tsv")

    def __call__(self):
        # run igblstn
        if self.seqtype == "TCR":
            self.igblast("TR")
        else:
            self.igblast("IG")
        
        self.get_metrics()

    @utils.add_log
    def igblast(self, chain):
        cwd = os.getcwd()
        # Igblastn program expects the internal_data directory under current directory
        os.chdir(self.soft_path)

        cmd = (
            f"./bin/igblastn "
            f"-query {self.fasta} "
            f"-organism {self.species} "
            f"-ig_seqtype {self.seqtype} "
            f"-auxiliary_data ./optional_file/{self.species}_gl.aux "
            f"-num_threads {self.args.thread} "
            f"-germline_db_V {self.ref_path}/{chain}V.fa "
            f"-germline_db_D {self.ref_path}/{chain}D.fa "
            f"-germline_db_J {self.ref_path}/{chain}J.fa "
            "-domain_system imgt -show_translation -outfmt 19 " # outfmt19 is an AIRR tab-delimited file, IgBLAST v1.9.0 or higher required.
            f"-out {self.airr_out} "
        )

        self.igblast.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
        os.chdir(cwd)

    @utils.add_log
    def get_metrics(self):
        df = pd.read_csv(self.airr_out, sep='\t')
        df.fillna("", inplace=True)
        total_reads = df.shape[0]
        
        df = df[(df["v_call"]!="") | ((df["d_call"]!="")) | ((df["j_call"]!=""))]
        self.add_metric(
            name="UMIs Mapped to Any VDJ Gene",
            value=df.shape[0],
            total=total_reads,
            help_info=f"UMI Mapped to any germline VDJ gene segments"
        )

        # UMIs with CDR3
        df_cdr3 = df[df["cdr3_aa"]!=""]
        self.add_metric(
            name=f"UMIs with CDR3",
            value=df_cdr3.shape[0],
            total=total_reads,
            help_info=f"UMIs with CDR3 sequence"
        )

        # UMIs with Correct CDR3
        df_correct_cdr3 = df_cdr3[~(df_cdr3["cdr3_aa"].str.contains(r"\*")) & ~(df_cdr3["cdr3_aa"].str.contains("X"))]
        self.add_metric(
            name="UMIs with CDR3",
            value=df_correct_cdr3.shape[0],
            total=total_reads,
            help_info=f"UMIs with CDR3 sequence"
        )
        
        # UMIs Mapped Confidently To VJ Gene
        df_confident = df[df["productive"]=="T"]
        self.add_metric(
            name=f"UMIs Mapped Confidently to VJ Gene",
            value=df_confident.shape[0],
            total=total_reads,
            help_info=f"UMIs mapped to VJ gene pairs and with correct CDR3"
        )

        # chains
        for chain in self.chains:
            df_chain = df_confident[df_confident.locus == chain]
            self.add_metric(
                name=f"UMIs Mapped to {chain}",
                value=df_chain.shape[0],
                total=total_reads,
                help_info=f"UMIs mapped confidently to {chain}"
            )


def mapping_vdj(args):
    with Mapping_vdj(args, display_title="Mapping") as runner:
        runner()


def get_opts_mapping_vdj(parser, sub_program):
    parser.add_argument('--soft_path', help='soft path for igblast')
    parser.add_argument('--ref_path', help='reference path for igblast')
    parser.add_argument("--type", help='TCR or BCR', required=True)
    parser.add_argument(
        '--species',
        choices=['human', 'mouse'],
        help='Default human. human or mouse.',
        default='human'
    )
    if sub_program:
        parser.add_argument(
            "--fasta",
            help="Required. Input fasta file.",
            required=True,
        )
        parser = s_common(parser)