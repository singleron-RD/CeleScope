import configparser
import subprocess

import celescope.tools.utils as utils
from celescope.tools.mkref import Mkref, parse_genomeDir
from celescope.tools.mkref import get_opts_mkref as opts


def parse_genomeDir_rna(genomeDir):
    return parse_genomeDir(genomeDir, entrys=('fasta', 'gtf', 'mt_gene_list'))


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

    def __init__(self, genome_type, args):
        Mkref.__init__(self, genome_type, args)
        self.fasta = args.fasta
        self.gtf = args.gtf
        self.mt_gene_list = args.mt_gene_list

        # out file
        self.refflat = f'{self.genome_name}.refFlat'

    @utils.add_log
    def build_star_index(self):
        cmd = (
            f'STAR \\\n'
            f'--runMode genomeGenerate \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--genomeDir {self.genomeDir} \\\n'
            f'--genomeFastaFiles {self.fasta} \\\n'
            f'--sjdbGTFfile {self.gtf} \\\n'
            f'--sjdbOverhang 100 \\\n'
        )
        Mkref_rna.build_star_index.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def write_config(self):
        config = configparser.ConfigParser()
        config['genome'] = {}
        genome = config['genome']
        genome['genome_name'] = self.genome_name
        genome['genome_type'] = self.genome_type
        genome['fasta'] = self.fasta
        genome['gtf'] = self.gtf
        genome['refFlat'] = self.refflat
        genome['mt_gene_list'] = self.mt_gene_list
        with open(self.config_file, 'w') as config_handle:
            config.write(config_handle)

    @utils.add_log
    def build_refflat(self):
        cmd = (
            'gtfToGenePred -genePredExt -geneNameAsName2 \\\n'
            f'{self.gtf} /dev/stdout | \\\n'
            'awk \'{print $12"\\t"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}\' \\\n'
            f'> {self.refflat} \\\n'
        )
        Mkref_rna.build_refflat.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        if not self.dry_run:
            self.build_refflat()
            self.build_star_index()
        self.write_config()


def mkref(args):
    genome_type = 'rna'
    runner = Mkref_rna(genome_type, args)
    runner.run()


def get_opts_mkref(parser, sub_program):
    opts(parser, sub_program)
    if sub_program:
        parser.add_argument(
            "--fasta",
            help="Required. Genome fasta file. Must be relative file path to genomeDir.",
            required=True
        )
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
