
from celescope.tools.mkref import Mkref, super_opts

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
    with Mkref_fusion(genome_type, args, files=('fusion_pos',), ) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    super_opts(parser, sub_program)
    if sub_program:
        parser.add_argument(
            "--fusion_pos",
            help="""
fusion position file. A two column tab-delimited text file with header.
"pos" is the end postion of the first gene(1-based).
e.g.  
name\tpos  
PML_3\t183  
PML_4\t254  
PML_5\t326  
PML_6\t204   
""",
            required=True,)

