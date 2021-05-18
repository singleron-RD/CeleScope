import celescope.tools.utils as utils
from celescope.tools.star import Star, get_opts_star


class StarVirus(Star):
    def __init__(self, args, step_name):
        Star.__init__(self, args, step_name)
        self.STAR_index = args.virus_genomeDir
        # change outPrefix
        self.outPrefix = f'{self.outdir}/{self.sample}_virus_'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'

    def run(self):
        self.STAR()
        self.sort_bam()
        self.index_bam()
        self.clean_up()


@utils.add_log
def star_virus(args):
    step_name = "star_virus"
    runner = StarVirus(args, step_name)
    runner.run()


def get_opts_star_virus(parser, sub_program):
    get_opts_star(parser, sub_program)
    parser.add_argument('--virus_genomeDir', help='virus genome dir', required=True)

