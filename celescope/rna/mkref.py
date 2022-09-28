from celescope.tools import utils
from celescope.tools.mkref import Mkref, super_opts


class Mkref_rna(Mkref):
    """
    ## Features
    - Create a genome reference directory.

    ## Usage
    ```
    celescope utils mkgtf Homo_sapiens.GRCh38.99.gtf Homo_sapiens.GRCh38.99.filtered.gtf
    celescope rna mkref \\
    --genome_name Homo_sapiens_ensembl_99_filtered \\
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \\
    --gtf Homo_sapiens.GRCh38.99.filtered.gtf
    ```

    ## Output

    - STAR genome index files

    - Genome config file
    ```
    $ cat celescope_genome.config
    [genome]
    genome_name = Homo_sapiens_ensembl_99
    genome_type = rna
    fasta = Homo_sapiens.GRCh38.dna.primary_assembly.fa
    gtf = Homo_sapiens.GRCh38.99.gtf
    ```
    """


    @utils.add_log
    def build_rna_star_index(self):
        cmd = (
            f'STAR \\\n'
            f'--runMode genomeGenerate \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--genomeDir ./ \\\n'
            f'--genomeFastaFiles {self.fasta} \\\n'
            f'--sjdbGTFfile {self.gtf} \\\n'
            f'--sjdbOverhang 100 \\\n'
        )
        if self.STAR_param:
            cmd += (" " + self.STAR_param)
        self.build_star_index.logger.info(cmd)
        self.debug_subprocess_call(cmd)


    @staticmethod
    def parse_genomeDir(genomeDir):
        return Mkref.parse_genomeDir(genomeDir, files=('gtf', 'mt_gene_list'))


    @utils.add_log
    def run(self):
        super().run()
        self.build_star_index()



def mkref(args):
    genome_type = 'rna'
    # files do not contain refflat because refflat is not input argument
    with Mkref_rna(genome_type, args, files=('gtf', 'mt_gene_list')) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    super_opts(parser, sub_program)
    if sub_program:
        parser.add_argument(
            "--gtf",
            help="Required. Genome gtf file. Use absolute path or relative path to `genomeDir`.",
            required=True
        )
        parser.add_argument(
            "--mt_gene_list",
            help="""Mitochondria gene list file. Use absolute path or relative path to `genomeDir`.
It is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.""",
            default="None"
        )
