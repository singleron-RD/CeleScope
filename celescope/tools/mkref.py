import abc
import configparser

from celescope.tools.__init__ import GENOME_CONFIG


def parse_genomeDir(genomeDir, entrys=None):
    config_file = f'{genomeDir}/{GENOME_CONFIG}'
    config = configparser.ConfigParser()
    config.read_file(open(config_file))
    genome = config['genome']
    # add path
    if entrys:
        for entry in entrys:
            if entry in genome and genome[entry] != 'None':
                genome[entry] = f'{genomeDir}/{genome[entry]}'
            else:
                genome[entry] = "None"
    return genome


class Mkref():

    def __init__(self, genome_type, args):
        self.genomeDir = args.genomeDir
        self.thread = args.thread
        self.genome_name = args.genome_name
        self.dry_run = args.dry_run
        self.genome_type = genome_type

        # out file
        self.config_file = f'{self.genomeDir}/{GENOME_CONFIG}'

    @abc.abstractmethod
    def run(self):
        return

    @abc.abstractmethod
    def write_config(self):
        return


def get_opts_mkref(parser, sub_program):
    if sub_program:
        parser.add_argument("--genomeDir", help="Default='./'. Output directory.", default='./')
        parser.add_argument("--thread", help="Default=6. Threads to use.", default=6)
        parser.add_argument("--genome_name", help="Required, genome name. ", required=True)
        parser.add_argument("--dry_run", help="Only write config file and exit.", action='store_true')
