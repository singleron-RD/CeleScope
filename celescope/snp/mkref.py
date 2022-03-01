import os
import subprocess

from celescope.tools import utils
from celescope.tools.mkref import Mkref, super_opts

class Mkref_snp(Mkref):
    """
    ## Features
    - Create dictionary file and fasta index for gatk SplitNCigarReads.
    (https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) 
    Need to run `celescope rna mkref` first

    ## Output
    - fasta index
    - gatk dictionary file

    ## Usage
    ```
    # run celescope rna mkref first
    celescope snp mkref \\
     --genome_name Homo_sapiens_ensembl_99 \\
     --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
    ```
    """


    @utils.add_log
    def build_fasta_index(self):
        cmd = (
            f'samtools faidx {self.fasta}'
        )
        self.build_fasta_index.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def build_fasta_dict(self):
        cmd = (
            f'gatk CreateSequenceDictionary -R {self.fasta}'
        )
        self.build_fasta_dict.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        assert os.path.exists(self.config_file), (
            f"{self.config_file} not exist! "
            "You need to build rna genome with 'celescope rna mkref' first."
        )
        super().run()
        self.build_fasta_index()
        self.build_fasta_dict()


def mkref(args):
    genome_type = 'snp'
    with Mkref_snp(genome_type, args, ) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    super_opts(parser, sub_program)
