from celescope.tools import utils
from celescope.tools.star_mixin import Star_mixin, get_opts_star_mixin
from celescope.capture_virus.mkref import Mkref_virus


class Star_virus(Star_mixin):
    """
    ## Features
    - Map reads to the viral genome using STAR.

    ## Output
    - `{sample}_virus_Aligned.sortedByCoord.out.bam` : Aligned BAM sorted by coordinate.
    """

    def __init__(self, args, display_title=None):
        # before init
        args.genomeDir = args.virus_genomeDir

        super().__init__(args, add_prefix='virus', display_title=display_title)
        self.genome = Mkref_virus.parse_genomeDir(self.genomeDir)


@utils.add_log
def star_virus(args):
    with Star_virus(args, display_title='Mapping') as runner:
        runner.run()


def get_opts_star_virus(parser, sub_program):
    get_opts_star_mixin(parser, sub_program)
    parser.add_argument('--virus_genomeDir', help='Virus genome dir.', required=True)
