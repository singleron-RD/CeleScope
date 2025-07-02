"""
bulk-seq vdj mapping
"""

import subprocess
import csv
import os

from celescope.tools import utils, step
from celescope.vdj import mapping_vdj as super_vdj
from celescope.vdj.__init__ import CHAINS
from celescope.bulk_rna.starsolo import get_barcode_sample
import celescope.tools.parse_chemistry as parse_chemistry
from celescope.__init__ import HELP_DICT
from celescope.chemistry_dict import chemistry_dict
import pandas as pd


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

        # out
        self.sample_productive_dir = f"{self.outdir}/sample_productive"
        if not os.path.exists(self.sample_productive_dir):
            os.makedirs(self.sample_productive_dir)
        self.airr_out = f"{self.out_prefix}_airr.tsv"
        self.productive_file = f"{self.out_prefix}_productive.tsv"

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
            barcode: open(f"{self.sample_productive_dir}/{sample}_productive.tsv", "wt")
            for barcode, sample in self.barcode_sample.items()
        }
        headers = [
            "sample",
            "barcode",
            "umi",
            "chain",
            "bestVGene",
            "bestDGene",
            "bestJGene",
            "nSeqCDR3",
            "aaSeqCDR3",
        ]
        for f in tsv_handles.values():
            f.write("\t".join(headers) + "\n")

        consensus_metrics_df = pd.read_csv(
            self.args.consensus_metrics_file, sep="\t", index_col=0
        )
        consensus_metrics_dict = consensus_metrics_df.to_dict(orient="index")
        metrics = utils.nested_defaultdict(dim=2)
        for barcode in consensus_metrics_dict:
            if barcode in self.barcode_sample:
                metrics[barcode]["sample"] = self.barcode_sample[barcode]
            else:
                metrics[barcode]["sample"] = ""
            metrics[barcode]["read"] = consensus_metrics_dict[barcode]["n_read"]

        total_umi = 0
        with open(self.airr_out, "rt") as infile:
            reader = csv.DictReader(infile, delimiter="\t")
            for row in reader:
                total_umi += 1
                barcode, umi = row["sequence_id"].split(":")[0:2]
                metrics[barcode]["UMI"] += 1
                if row["v_call"] != "" or row["d_call"] != "" or row["j_call"] != "":
                    metrics[barcode]["UMI_mapped_to_any_vdj"] += 1
                else:
                    continue
                if row["productive"] == "T" and "N" not in row["junction"]:
                    metrics[barcode]["UMI_confident"] += 1
                    metrics[barcode][f"UMI_confident_{row['locus']}"] += 1
                else:
                    continue
                if barcode in self.barcode_sample:
                    line = [
                        self.barcode_sample[barcode],
                        barcode,
                        umi,
                        row["locus"],
                        row["v_call"],
                        row["d_call"],
                        row["j_call"],
                        row["junction"],
                        row["junction_aa"],
                    ]
                    tsv_handles[barcode].write("\t".join(line) + "\n")

        umi_mapped_to_any_vdj = sum(
            v for barcode, k in metrics.items() for v in [k["UMI_mapped_to_any_vdj"]]
        )
        self.add_metric(
            name="UMIs Mapped to Any VDJ Gene",
            value=umi_mapped_to_any_vdj,
            total=total_umi,
            help_info="UMIs Mapped to any germline VDJ gene segments",
        )

        umi_confident = sum(
            v for barcode, k in metrics.items() for v in [k["UMI_confident"]]
        )
        self.add_metric(
            name="UMIs Mapped Confidently to VJ Gene",
            value=umi_confident,
            total=total_umi,
            help_info="UMIs with productive rearrangement mapped to VJ gene pairs and without N in the junction sequence",
        )

        for chain in self.chains:
            UMI_confident_chain = sum(
                v
                for barcode, k in metrics.items()
                for v in [k[f"UMI_confident_{chain}"]]
            )
            self.add_metric(
                name=f"UMIs Mapped Confidently to {chain}",
                value=UMI_confident_chain,
                total=total_umi,
                help_info=f"UMIs with productive rearrangement mapped to {chain}",
            )

        metrics_df = pd.DataFrame.from_dict(metrics, orient="index")
        metrics_df.index.name = "barcode"
        metrics_df.reset_index(inplace=True)
        metrics_df = metrics_df.sort_values(by="read", ascending=False)

        metrics_df.to_csv(f"{self.out_prefix}_well_metrics.tsv", sep="\t", index=False)

    def run(self):
        # run igblstn
        # self.igblast()
        self.process_airr()


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
