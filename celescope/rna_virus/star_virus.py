import celescope.tools.utils as utils
from celescope.tools.star_mixin import StarMixin, get_opts_star_mixin
from celescope.tools.step import Step


class StarVirus(Step, StarMixin):
    """
    star virus class
    """
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        # add genomeDir 
        args.genomeDir = args.virus_genomeDir
        
        StarMixin.__init__(self, args)

        # change outPrefix
        self.outPrefix = f'{self.outdir}/{self.sample}_virus_'
        self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
        self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'

    def run(self):
        self.run_star()
        self.clean_up()


@utils.add_log
def star_virus(args):
    step_name = "star_virus"
    runner = StarVirus(args, step_name)
    runner.run()


def get_opts_star_virus(parser, sub_program):
    get_opts_star_mixin(parser, sub_program)
    parser.add_argument('--virus_genomeDir', help='virus genome dir', required=True)

