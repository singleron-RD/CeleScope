import celescope.tools.utils as utils
from celescope.tools.star import Star, get_opts_star


class StarFusion(Star):
    def __init__(self, args, step_name):
        Star.__init__(self, args, step_name)
        self.STAR_index = args.fusion_genomeDir

    def run(self):
        self.STAR()
        self.sort_bam()
        self.index_bam()
        self.clean_up()


@utils.add_log
def star_fusion(args):
    step_name = "star_fusion"
    runner = StarFusion(args, step_name)
    runner.run()


def get_opts_star_fusion(parser, sub_program):
    get_opts_star(parser, sub_program)
    # will cause `conflicting option string: --genomeDir`
    # parser.add_argument('--genomeDir', help=argparse.SUPPRESS)
    parser.add_argument('--fusion_genomeDir', help='fusion gene STAR index genome directory', required=True)

