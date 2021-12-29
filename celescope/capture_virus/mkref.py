
from celescope.tools.mkref import Mkref, super_opts

class Mkref_virus(Mkref):
    """
    Features
    - Create a virus genome reference directory.

    Output

    - STAR genome index files
    - Genome config file

    Usage
    ```
    celescope capture_virus mkref \\
        --genome_name EBV \\
        --fasta EBV_genome.fasta \\
        --genomeSAindexNbases 7
    ```

    ```
    $ cat celescope_genome.config
    [genome]
    genome_type = virus
    fasta = EBV_genome.fasta
    genome_name = EBV
    genomesaindexnbases = 7
    ```
    """

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
