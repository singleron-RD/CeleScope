import configparser
import subprocess

import celescope.tools.utils as utils
from celescope.tools.mkref import Mkref, super_opts

class Mkref_virus(Mkref):

    def run(self):
        super().run()
        self.build_star_index()

def mkref(args):
    genome_type = 'virus'
    with Mkref_virus(genome_type, args, non_files=('genomeSAindexNbases',)) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    super_opts(parser, sub_program)
    parser.add_argument("--genomeSAindexNbases", help="STAR genomeSAindexNbases", default=14)
