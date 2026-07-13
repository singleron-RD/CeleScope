import subprocess
import sys

from pathlib import Path

from celescope.tools.starsolo import (
    Starsolo as tools_Starsolo,
    get_opts_starsolo as tools_opts,
)
from celescope.tools.matrix import CountMatrix
from celescope.tools import utils


class StarsoloMirna(tools_Starsolo):
    """
    Starsolo for miRNA alignment
    """

    def __init__(self, args):
        # Call parent init but override genomeDir with mirna_genomeDir
        # Save original genomeDir
        self.original_genomeDir = args.genomeDir
        args.genomeDir = args.mirna_genomeDir
        super().__init__(args)

        # Add miRNA-specific STAR parameters
        self.extra_starsolo_args += (
            " --outFilterMismatchNoverLmax 0.1"
            " --outFilterMatchNmin 15"
            " --outFilterScoreMinOverLread 0.9"
            " --outFilterMatchNminOverLread 0.9"
            " --alignIntronMax 1"
        )

        self.mrna_raw = CountMatrix.from_matrix_dir(f"{self.outs_dir}/raw/")
        self.mrna_filtered = CountMatrix.from_matrix_dir(f"{self.outs_dir}/filtered/")

        # output files
        self.solo_out_dir = f"{self.outdir}/{self.sample}_Solo.out/"
        solo_dir = f"{self.outdir}/{self.sample}_Solo.out/{args.report_soloFeature}"
        self.raw_matrix = Path(f"{solo_dir}/raw")

        self.outs = []

    @utils.add_log
    def concat_matrix(self):
        mirna_raw = CountMatrix.from_matrix_dir(self.raw_matrix)
        concat_raw = self.mrna_raw.concat_by_barcodes(mirna_raw)
        concat_raw.to_matrix_dir(f"{self.outs_dir}/raw/")

        concat_filtered = concat_raw.slice_matrix_bc(self.mrna_filtered.get_barcodes())
        concat_filtered.to_matrix_dir(f"{self.outs_dir}/filtered/")

    @utils.add_log
    def gzip_matrix(self):
        cmd = f"gzip {self.raw_matrix}/*"
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_starsolo()
        self.gzip_matrix()
        self.concat_matrix()


def mirna(args):
    if not args.mirna_genomeDir:
        sys.stderr.write("Skip miRNA analysis.\n")
        return

    with StarsoloMirna(args) as runner:
        runner.run()


def get_opts_mirna(parser, sub_program=True):
    tools_opts(parser, sub_program)
    parser.add_argument(
        "--mirna_genomeDir",
        help="miRNA genome directory for STAR alignment",
    )
    parser.set_defaults(
        outFilterMatchNmin=15,
        soloCellFilter="None",
    )
    return parser
