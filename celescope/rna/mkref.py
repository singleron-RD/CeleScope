import subprocess

import celescope.tools.utils as utils
from celescope.tools.mkref import Mkref, super_opts


class Mkref_rna(Mkref):
    """
    Features
    - Create a genome reference directory.

    Output

    - STAR genome index files

    - Genome refFlat file

    - Genome config file
    ```
    $ cat celescope_genome.config
    [genome]
    genome_name = Homo_sapiens_ensembl_99
    genome_type = rna
    fasta = Homo_sapiens.GRCh38.dna.primary_assembly.fa
    gtf = Homo_sapiens.GRCh38.99.gtf
    refflat = Homo_sapiens_ensembl_99.refFlat
    ```
    """

    @utils.add_log
    def build_star_index(self):
        cmd = (
            f'STAR \\\n'
            f'--runMode genomeGenerate \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--genomeDir ./ \\\n'
            f'--genomeFastaFiles {self.fasta} \\\n'
            f'--sjdbGTFfile {self.gtf} \\\n'
            f'--sjdbOverhang 100 \\\n'
        )
        self.build_star_index.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


    @utils.add_log
    def build_refflat(self):
        cmd = (
            'gtfToGenePred -genePredExt -geneNameAsName2 \\\n'
            f'{self.gtf} /dev/stdout | \\\n'
            'awk \'{print $12"\\t"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}\' \\\n'
            f'> {self.refflat} \\\n'
        )
        self.build_refflat.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    def parse_genomeDir(genomeDir):
        super().parse_genomeDir(genomeDir, files=('gtf', 'mt_gene_list'))

    @utils.add_log
    def run(self):
        super().run()
        self.build_star_index()
        self.build_refflat()


def mkref(args):
    genome_type = 'rna'
    with Mkref_rna(genome_type, args, files=('gtf', 'mt_gene_list'), non_files=('genomeSAindexNbases',)) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    super_opts(parser, sub_program)
    if sub_program:
        parser.add_argument(
            "--gtf",
            help="Required. Genome gtf file. Must be relative file path to genomeDir.",
            required=True
        )
        parser.add_argument(
            "--mt_gene_list",
            help="""Mitochondria gene list file. Must be relative file path to genomeDir.
It is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.""",
            default="None"
        )
        parser.add_argument("--genomeSAindexNbases", help="STAR genomeSAindexNbases", default=14)
