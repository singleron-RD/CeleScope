import subprocess
import sys

from celescope.tools import utils
from celescope.tools.make_ref import MakeRef_STAR


def parse_attributes(attrs):
    dic = {}
    for attr_str in attrs.split(";"):
        if attr_str:
            attr, val = attr_str.split("=")
            val = set(val.split(","))
            dic[attr] = val
    return dic


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
        parser.add_argument(
            "--attributes",
            help="Attributes to keep. Example: `gene_biotype=protein_coding,lncRNA,antisense;`",
            default="gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene;",
        )
