import configparser
import os
import subprocess

import celescope.tools.utils as utils
from celescope.tools.mkref import Mkref
from celescope.tools.mkref import get_opts_mkref as opts
from celescope.__init__ import HELP_DICT


class Mkref_snp(Mkref):
    """
    Features
    - Create dictionary file and fasta index for gatk SplitNCigarReads.
    (https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) 
    Need to run `celescope rna mkref` first

    Output
    - fasta index
    - gatk dictionary file

    Usage
    ```
    # run celescope rna mkref first
    celescope snp mkref \\
     --genome_name Homo_sapiens_ensembl_99 \\
     --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
    ```
    """

    def __init__(self, genome_type, args):
        Mkref.__init__(self, genome_type, args)
        self.fasta = args.fasta

    @utils.add_log
    def write_config(self):
        config = configparser.ConfigParser()
        config.read(self.config_file)
        genome = config['genome']
        genome['fasta_gatk_index'] = str(True)
        with open(self.config_file, 'w') as config_handle:
            config.write(config_handle)

    @utils.add_log
    def build_fasta_index(self):
        cmd = (
            f'samtools faidx {self.fasta}'
        )
        Mkref_snp.build_fasta_index.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def build_fasta_dict(self):
        cmd = (
            f'gatk CreateSequenceDictionary -R {self.fasta}'
        )
        Mkref_snp.build_fasta_dict.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        os.chdir(self.genomeDir)
        assert os.path.exists(self.config_file), (
            f"{self.config_file} not exist! "
            "You need to build rna genome with 'celescope rna mkref' first."
        )
        if not self.dry_run:
            self.build_fasta_index()
            self.build_fasta_dict()
        self.write_config()


def mkref(args):
    genome_type = 'snp'
    runner = Mkref_snp(genome_type, args)
    runner.run()


def get_opts_mkref(parser, sub_program):
    opts(parser, sub_program)
    if sub_program:
        parser.add_argument("--fasta", help=HELP_DICT['fasta'], required=True)
