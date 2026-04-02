import subprocess
import sys

from celescope.tools import utils
from celescope.tools.make_ref import MakeRef_STAR


class Mkref_rna(MakeRef_STAR):
    def __init__(self, genome_type, args):
        super().__init__(genome_type, args)
        self.files["gtf"] = self.args.gtf
        self.files["mt_gene_list"] = self.args.mt_gene_list

    @utils.add_log
    def build_rna_star_index(self):
        SA = self._get_SA()
        cmd = (
            f"STAR \\\n"
            f"--runMode genomeGenerate \\\n"
            f"--runThreadN {self.args.thread} \\\n"
            f"--genomeDir ./ \\\n"
            f"--genomeFastaFiles {self.args.fasta} \\\n"
            f"--sjdbGTFfile {self.args.gtf} \\\n"
            f"--sjdbOverhang 100 \\\n"
            f"--genomeSAindexNbases {SA} \\\n"
        )
        if self.STAR_param:
            cmd += " " + self.args.STAR_param
        sys.stderr.write(cmd + "\n")
        if not self.dry_run:
            subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        self.build_rna_star_index()


def mkref(args):
    genome_type = "rna"
    with Mkref_rna(genome_type, args) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    MakeRef_STAR.opts(parser, sub_program)
    if sub_program:
        parser.add_argument("--gtf", help="Required. Gtf file name.", required=True)
        parser.add_argument(
            "--mt_gene_list",
            help="""Mitochondria gene list file name. This file is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.""",
            default="None",
        )
