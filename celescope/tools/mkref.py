import configparser
import subprocess
import sys

import celescope.tools.utils as utils
from celescope.tools.__init__ import GENOME_CONFIG
from celescope.__init__ import HELP_DICT




class Mkref():

    def __init__(self, genome_type, args, files=(), non_files=()):
        self.thread = args.thread
        self.genome_type = genome_type
        self.dry_run = args.dry_run
        self.files = ('fasta',) + files
        self.non_files = ('genome_name',) + non_files

        for entry in self.files + self.non_files:
            value = str(getattr(args, entry))
            setattr(self, entry, value)

        self.config_dict = {}

        # out file
        self.config_file = GENOME_CONFIG

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
        self.build_star_index.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def set_config_dict(self):
        """
        if some files are not in input arguments, set them in overwrite set_config_dict function
        """
        self.config_dict['genome_type'] = self.genome_type

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

        >>> raw_file_path = '/root/Homo_sapiens.GRCh38.92.chr.gtf'
        >>> Mkref.get_file_path(raw_file_path, 'fake')
        '/root/Homo_sapiens.GRCh38.92.chr.gtf'
        """
        file_path = raw_file_path
        if file_path and not file_path.startswith('/') and file_path != "None":
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
        genome = config['genome']

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
