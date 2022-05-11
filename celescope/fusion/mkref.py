
from celescope.tools.mkref import Mkref, super_opts
from celescope.__init__ import HELP_DICT

class Mkref_fusion(Mkref):
    """
    ## Features
    - Create a fusion genome directory.

    ## Output

    - STAR genome index files
    - Genome config file

    ## Usage
    ```
    celescope fusion mkref \\
    --genome_name {genome_name} \\
    --fasta fusion.fasta \\
    --fusion_pos fusion_pos.txt \\
    --genomeSAindexNbases 4
    ```
    """

    def run(self):
        super().run()
        self.build_star_index()

    @staticmethod
    def parse_genomeDir(genomeDir):
        return Mkref.parse_genomeDir(genomeDir, files=('fusion_pos',))



def mkref(args):
    genome_type = 'fusion'
    with Mkref_fusion(genome_type, args, files=('fusion_pos',), non_files=('genomeSAindexNbases',)) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    super_opts(parser, sub_program)
    if sub_program:
        parser.add_argument("--fusion_pos",help=HELP_DICT['fusion_pos'],required=True,)
        parser.add_argument("--genomeSAindexNbases", help=HELP_DICT['genomeSAindexNbases'], default=14)
