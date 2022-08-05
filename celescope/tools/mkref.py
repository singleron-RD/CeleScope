import configparser
import subprocess
import sys
import os
import math

from celescope.tools import utils
from celescope.tools.__init__ import GENOME_CONFIG
from celescope.__init__ import HELP_DICT, __VERSION__


class Mkref():

    def __init__(self, genome_type, args, files=(), non_files=()):
        self.thread = args.thread
        self.genome_type = genome_type
        self.celescope_version = __VERSION__
        self.dry_run = args.dry_run
        self.STAR_param = args.STAR_param
        self.files = ('fasta',) + files
        self.non_files = ('genome_name',) + non_files

        for entry in self.files + self.non_files:
            value = str(getattr(args, entry))
            setattr(self, entry, value)

        self.config_dict = {}

        # get genomeSAindexNbases from fasta size
        fasta_size = os.path.getsize(self.fasta)
        self.genomeSAindexNbases = Mkref.get_SA(fasta_size)
        print(f'genomeSAindexNbases: {self.genomeSAindexNbases}')

        # out file
        self.config_file = GENOME_CONFIG

    @staticmethod
    def get_SA(fasta_size):
        """
        Get genomeSAindexNbases from fasta size. For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical value of min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 8, for 100 kiloBase genome, this is equal to 7.

        As mentioned in STAR manual, or 1 megaBase genome, this is equal to 9. But when the genome is relatively small, e.g. < 3kb,
        using large value may caused core dump error when mapping.

        Args:
            fasta_size: fasta size in bytes

        >>> Mkref.get_SA(10 ** 6)
        8
        >>> Mkref.get_SA(10 ** 5)
        7
        >>> Mkref.get_SA(10 ** 50)
        14
        """
        return min(int(math.log2(fasta_size) / 2 - 1), 14)

    def run(self):
        if self.dry_run:
            sys.exit(0)

    @utils.add_log
    def build_star_index(self):
        cmd = (
            f'STAR \\\n'
            f'--runMode genomeGenerate \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--genomeDir ./ \\\n'
            f'--genomeFastaFiles {self.fasta} \\\n'
            f'--genomeSAindexNbases {self.genomeSAindexNbases} \\\n'
        )
        if self.STAR_param:
            cmd += (" " + self.STAR_param)
        self.build_star_index.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def set_config_dict(self):
        """
        if some files are not in input arguments, set them in overwrite set_config_dict function
        """
        self.config_dict['genome_type'] = self.genome_type
        self.config_dict['celescope_version'] = self.celescope_version

        for entry in self.files + self.non_files:
            value = getattr(self, entry)
            print(f'{entry} : {value}')
            self.config_dict[entry] = value

    @utils.add_log
    def _write_config(self):
        config = configparser.ConfigParser()
        config['genome'] = self.config_dict

        with open(self.config_file, 'w') as config_handle:
            config.write(config_handle)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.set_config_dict()
        self._write_config()

    @staticmethod
    def get_file_path(raw_file_path, genomeDir):
        """
        if raw_file_path is not absolute path and not None,
        add genomeDir
        if 'NONE'(str), return None

        >>> raw_file_path = '/root/Homo_sapiens.GRCh38.92.chr.gtf'
        >>> Mkref.get_file_path(raw_file_path, 'fake')
        '/root/Homo_sapiens.GRCh38.92.chr.gtf'

        >>> raw_file_path = 'None'
        >>> Mkref.get_file_path(raw_file_path, 'fake')
        """
        if (not raw_file_path) or (raw_file_path == 'None'):
            return None

        file_path = raw_file_path
        if not file_path.startswith('/'):
            file_path = f'{genomeDir}/{file_path}'

        return file_path            

    @staticmethod
    def parse_genomeDir(genomeDir, files=()):
        """
        Parse genomeDir and return a dict.
        """
        files = ('fasta',) + files

        config_file = f'{genomeDir}/{GENOME_CONFIG}'
        config = configparser.ConfigParser()
        config.read_file(open(config_file))
        genome = dict(config['genome'])

        for entry in files:
            if entry not in genome:
                raise ValueError(f'{entry} not in {config_file}')
            genome[entry] = Mkref.get_file_path(genome[entry], genomeDir)        

        return genome


def super_opts(parser, sub_program):
    if sub_program:
        parser.add_argument("--thread", help="Default=6. Threads to use.", default=6)
        parser.add_argument("--genome_name", help="Required, genome name. ", required=True)
        parser.add_argument("--dry_run", help="Only write config file and exit.", action='store_true')
        parser.add_argument("--fasta", help=HELP_DICT['fasta'], required=True)
        parser.add_argument('--STAR_param', help=HELP_DICT['additional_param'], default="")
