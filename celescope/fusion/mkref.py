import configparser
import subprocess

import celescope.tools.utils as utils
from celescope.tools.mkref import Mkref
from celescope.tools.mkref import get_opts_mkref as opts
from celescope.tools.mkref import parse_genomeDir


def parse_genomeDir_fusion(genomeDir):
    return parse_genomeDir(genomeDir, entrys=('fasta', 'fusion_pos'))


class Mkref_fusion(Mkref):
    """
    Features
    - Create a fusion genome directory.

    Output

    - STAR genome index files
    - Genome config file
    
    Usage
    ```
    celescope fusion mkref \\
    --genome_name {genome_name} \\
    --fasta fusion.fasta \\
    --fusion_pos fusion_pos.txt \\
    --genomeSAindexNbases 4
    ```
    """
    def __init__(self, genome_type, args):
        Mkref.__init__(self, genome_type, args)
        self.fasta = args.fasta
        self.genomeSAindexNbases = args.genomeSAindexNbases
        self.fusion_pos = args.fusion_pos

    @utils.add_log
    def build_star_index(self):
        cmd = (
            f'STAR \\\n'
            f'--runMode genomeGenerate \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--genomeDir {self.genomeDir} \\\n'
            f'--genomeFastaFiles {self.fasta} \\\n'
            f'--genomeSAindexNbases {self.genomeSAindexNbases} \\\n'
        )
        Mkref_fusion.build_star_index.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def write_config(self):
        config = configparser.ConfigParser()
        config['genome'] = {}
        genome = config['genome']
        genome['genome_name'] = self.genome_name
        genome['genome_type'] = self.genome_type
        genome['fasta'] = self.fasta
        genome['fusion_pos'] = self.fusion_pos
        genome['genomeSAindexNbases'] = self.genomeSAindexNbases
        with open(self.config_file, 'w') as config_handle:
            config.write(config_handle)

    def run(self):
        if not self.dry_run:
            self.build_star_index()
        self.write_config()


def mkref(args):
    genome_type = 'fusion'
    runner = Mkref_fusion(genome_type, args)
    runner.run()


def get_opts_mkref(parser, sub_program):
    opts(parser, sub_program)
    if sub_program:
        parser.add_argument("--fasta", help="Fusion fasta file.", required=True)
        parser.add_argument(
            "--fusion_pos",
            help="""
fusion position file. A two column tab-delimited text file with header.
"pos" is the end postion of the first gene(1-based).
e.g.  
tag\tpos  
PML_3\t183  
PML_4\t254  
PML_5\t326  
PML_6\t204   
""",
            required=True,)
        parser.add_argument("--genomeSAindexNbases", help="STAR genomeSAindexNbases", default=4)
