import os
import subprocess
import sys

from celescope.tools import utils
from celescope.tools.make_ref import MakeRef, GENOME_CONFIG


class Mkref_snp:
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

    def __init__(self):
        pass

    @utils.add_log
    def build_fasta_index(self, fasta):
        cmd = f"samtools faidx {fasta}"
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def build_fasta_dict(self, fasta):
        cmd = f"gatk CreateSequenceDictionary -R {fasta}"
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        assert os.path.exists(GENOME_CONFIG), (
            f"{GENOME_CONFIG} not exist! "
            "You need to build rna genome with 'celescope rna mkref' first."
        )
        # reset genome type here
        config = MakeRef.get_config("./")
        fasta = config["files"]["fasta"]
        self.build_fasta_index(fasta)
        self.build_fasta_dict(fasta)

        config["meta"]["genome_type"] = "rna, snp"
        with open(GENOME_CONFIG, "w") as config_handle:
            config.write(config_handle)


def mkref(args):
    # pylint: disable=unused-argument
    Mkref_snp().run()


def get_opts_mkref(parser, sub_program):
    # pylint: disable=unused-argument
    pass
