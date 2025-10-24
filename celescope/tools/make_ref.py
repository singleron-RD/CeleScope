import configparser
import subprocess
import sys
import os
import math

from celescope.tools.__init__ import GENOME_CONFIG
from celescope.__init__ import HELP_DICT, __VERSION__
from celescope.tools import utils


class MakeRef:
    def __init__(self, genome_type, args):
        self.config_file = GENOME_CONFIG

        self.args = args
        self.thread = args.thread
        self.dry_run = args.dry_run

        self.files = {}
        self.meta = {}
        self.meta["genome_name"] = args.genome_name
        self.meta["genome_type"] = genome_type
        self.meta["celescope_version"] = __VERSION__

    def write_config(self):
        config = configparser.ConfigParser()
        config.optionxform = str
        config["files"] = self.files
        config["meta"] = self.meta

        with open(GENOME_CONFIG, "w") as config_handle:
            config.write(config_handle)

    @staticmethod
    def get_config(genomeDir):
        """
        add genomeDir prefix to files
        """
        config_file = f"{genomeDir}/{GENOME_CONFIG}"
        if not os.path.exists(config_file):
            sys.exit(
                f"Error: {config_file} not found.\n"
                "Solution: Use the mkref command of CeleScope to generate the genome.\n"
            )
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(config_file)
        if "meta" not in config:
            sys.exit(
                "Error: CeleScope version >= 2.0.0 have upgraded the STAR version, so the old genome can no longer be used.\n"
                "Solution: Use the mkref command of CeleScope to regenerate the genome.\n"
                f"Genome path: {genomeDir}\n"
            )
        for k, v in dict(config["files"]).items():
            if v and v != "None":
                config["files"][k] = f"{genomeDir}/{v}"
            else:
                config["files"][k] = ""
        return config

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.write_config()

    @staticmethod
    def opts(parser, sub_program):
        if sub_program:
            parser.add_argument(
                "--thread", help="Default=6. Threads to use.", default=16
            )
            parser.add_argument(
                "--genome_name", help="Required, genome name. ", required=True
            )
            parser.add_argument(
                "--dry_run",
                help="Only write config file and exit.",
                action="store_true",
            )


class MakeRef_STAR(MakeRef):
    def __init__(self, genome_type, args):
        super().__init__(genome_type, args)
        self.fasta = args.fasta
        self.STAR_param = args.STAR_param
        self.files["fasta"] = args.fasta
        self.meta["STAR_param"] = args.STAR_param

    @staticmethod
    def opts(parser, sub_program):
        MakeRef.opts(parser, sub_program)
        if sub_program:
            parser.add_argument(
                "--fasta", help="Required. fasta file name.", required=True
            )
            parser.add_argument(
                "--STAR_param", help=HELP_DICT["additional_param"], default=""
            )

    @staticmethod
    def get_SA(fasta_size):
        """
        Get genomeSAindexNbases from fasta size. For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 8, for 100 kiloBase genome, this is equal to 7.

        As mentioned in STAR manual, or 1 megaBase genome, this is equal to 9. But when the genome is relatively small, e.g. < 3kb,
        using large value may caused core dump error when mapping.

        Args:
            fasta_size: fasta size in bytes

        >>> MakeRef_STAR.get_SA(10 ** 6)
        8
        >>> MakeRef_STAR.get_SA(10 ** 5)
        7
        >>> MakeRef_STAR.get_SA(10 ** 50)
        14
        """
        return min(int(math.log2(fasta_size) / 2 - 1), 14)

    def _get_SA(self):
        fasta_size = os.path.getsize(self.fasta)
        SA = self.get_SA(fasta_size)
        self.meta["genomeSAindexNbases"] = SA
        return SA

    @utils.add_log
    def build_star_index(self):
        SA = self._get_SA()
        cmd = (
            f"STAR \\\n"
            f"--runMode genomeGenerate \\\n"
            f"--runThreadN {self.thread} \\\n"
            f"--genomeDir ./ \\\n"
            f"--genomeFastaFiles {self.fasta} \\\n"
            f"--genomeSAindexNbases {SA} \\\n"
        )
        if self.STAR_param:
            cmd += " " + self.args.STAR_param
        sys.stderr.write(cmd + "\n")
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        self.build_star_index()
