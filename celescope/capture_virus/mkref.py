
from celescope.tools.mkref import Mkref, super_opts

class Mkref_virus(Mkref):
    """
    ## Features
    - Create a virus genome reference directory.

    ## Output

    - STAR genome index files
    - Genome config file

    ## Usage
    ```
    celescope capture_virus mkref \\
        --genome_name EBV \\
        --fasta EBV_genome.fasta \\
    ```

    ```
    $ cat celescope_genome.config
    [genome]
    genome_type = virus
    fasta = EBV_genome.fasta
    genome_name = EBV
    ```
    """

    def run(self):
        super().run()
        self.build_star_index()

def mkref(args):
    genome_type = 'virus'
    with Mkref_virus(genome_type, args, ) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    super_opts(parser, sub_program)
