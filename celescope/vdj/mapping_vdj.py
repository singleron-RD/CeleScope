'''
vdj mapping
'''
import pandas as pd
import subprocess

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj.__init__ import CHAINS


class Mapping_vdj(Step):
    """
    ## Features
    - Align R2 reads to IGMT(http://www.imgt.org/) database sequences with blast.

    ## Output
    - `{sample}_airr.tsv` The alignment result of each UMI.
    A tab-delimited file compliant with the AIRR Rearrangement schema(https://docs.airr-community.org/en/stable/datarep/rearrangements.html)

    - `{sample}_UMI_count_unfiltered.tsv` UMI reading for each (barcode, chain, VJ_pair) combination.

    - `{sample}_UMI_count_filtered.tsv` For each (barcode, chain) combination, only the record with the 
    most VJ_pair UMI reads is kept.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.ref_path = args.ref_path
        self.seqtype = args.type
        self.species = args.species
        self.chains = CHAINS[self.seqtype]

        # out
        self.airr_out = f'{self.out_prefix}_airr.tsv'
        self.UMI_count_unfiltered_file = f'{self.out_prefix}_UMI_count_unfiltered.tsv'
        self.UMI_count_filtered_file = f'{self.out_prefix}_UMI_count_filtered.tsv'

    def run(self):
        # run igblstn
        if self.seqtype == "TCR":
            self.igblast(chain="TR", ig_seqtype="TCR")
        else:
            self.igblast(chain="IG", ig_seqtype="Ig")
        self.mapping_summary()


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
            "-domain_system imgt -show_translation -outfmt 19 " # outfmt19 is an AIRR tab-delimited file, IgBLAST v1.9.0 or higher required.
            f"-out {self.airr_out} "
        )

        self.igblast.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def mapping_summary(self):
        df = pd.read_csv(self.airr_out, sep='\t')
        df.fillna("", inplace=True)
        total_reads = df.shape[0]

        self.add_metric(
            name="Species",
            value=self.species,
            help_info="Human or Mouse"
        )
        
        df = df[(df["v_call"]!="") | ((df["d_call"]!="")) | ((df["j_call"]!=""))]
        self.add_metric(
            name="UMIs Mapped to Any VDJ Gene",
            value=df.shape[0],
            total=total_reads,
            help_info="UMI Mapped to any germline VDJ gene segments"
        )

        # UMIs with CDR3
        df_cdr3 = df[df["cdr3_aa"]!=""]
        self.add_metric(
            name="UMIs with CDR3",
            value=df_cdr3.shape[0],
            total=total_reads,
            help_info="UMIs with CDR3 sequence"
        )

        # UMIs with Correct CDR3
        df_correct_cdr3 = df_cdr3[~(df_cdr3["cdr3_aa"].str.contains(r"\*")) & ~(df_cdr3["cdr3_aa"].str.contains("X"))]
        self.add_metric(
            name="UMIs with Correct CDR3",
            value=df_correct_cdr3.shape[0],
            total=total_reads,
            help_info="UMIs with CDR3 might have stop codon and these UMIs(or Reads) are classified as incorrect"
        )
        
        # UMIs Mapped Confidently To VJ Gene
        df_confident = df_correct_cdr3[df_correct_cdr3["productive"]=="T"]
        self.add_metric(
            name="UMIs Mapped Confidently to VJ Gene",
            value=df_confident.shape[0],
            total=total_reads,
            help_info="UMIs with productive rearrangement mapped to VJ gene pairs and with correct CDR3"
        )

        # UMIs Mapped Confidently to each chain
        for chain in self.chains:
            df_chain = df_confident[df_confident.locus == chain]
            self.add_metric(
                name=f"UMIs Mapped to {chain}",
                value=df_chain.shape[0],
                total=total_reads,
                help_info=f"UMIs mapped confidently to {chain}"
            )

        # output file
        df_VJ = df_confident[["sequence_id","locus", "v_call", "d_call", "j_call", "junction", "junction_aa"]]
        df_VJ.rename(columns={
            "sequence_id": "readID",
            "locus": "chain",
            "v_call": "bestVGene",
            "d_call": "bestDGene",
            "j_call": "bestJGene",
            "junction": "nSeqCDR3",
            "junction_aa": "aaSeqCDR3"
        }, inplace=True)
        df_VJ["barcode"] = df_VJ["readID"].apply(lambda x: x.split("_")[0])
        df_VJ["UMI"] = df_VJ["readID"].apply(lambda x: x.split("_")[1])

        for i in ["bestVGene", "bestDGene", "bestJGene"]:
            df_VJ[i] = df_VJ[i].apply(lambda x: x.split("*")[0])

        # filter1: keep top 1 in each combinations
        groupby_elements = [
            'barcode',
            'chain',
            'bestVGene',
            'bestJGene',
            'aaSeqCDR3',
            'nSeqCDR3',
        ]

        # unique UMI
        df_UMI = df_VJ.drop_duplicates(subset=["barcode", "UMI"], keep="first")
        df_UMI_count = df_UMI.groupby(
            groupby_elements, as_index=False).agg({"UMI": "count"})
        df_UMI_count = df_UMI_count.sort_values("UMI", ascending=False)
        # out unfiltered
        df_UMI_count.to_csv(self.UMI_count_unfiltered_file, sep="\t", index=False)

        df_UMI_count_filter = df_UMI_count.groupby(
            ["barcode", "chain"], as_index=False).head(1)
        # out filtered
        df_UMI_count_filter.to_csv(
            self.UMI_count_filtered_file,
            sep="\t",
            index=False
        )


def mapping_vdj(args):
    with Mapping_vdj(args, display_title="Mapping") as runner:
        runner.run()


def get_opts_mapping_vdj(parser, sub_program):
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