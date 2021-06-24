import celescope.tools.utils as utils
from celescope.tools.star_mixin import StarMixin, get_opts_star_mixin
from celescope.tools.step import Step


class StarFusion(Step, StarMixin):
    def __init__(self, args, step_name):
        args.genomeDir = args.fusion_genomeDir
        Step.__init__(self, args, step_name)
        StarMixin.__init__(self, args)

    def run(self):
        self.run_star()
        self.clean_up()


@utils.add_log
def star_fusion(args):
    step_name = "star_fusion"
    runner = StarFusion(args, step_name)
    runner.run()


def get_opts_star_fusion(parser, sub_program):
    get_opts_star_mixin(parser, sub_program)
    # will cause `conflicting option string: --genomeDir`
    # parser.add_argument('--genomeDir', help=argparse.SUPPRESS)
    parser.add_argument('--fusion_genomeDir', help='fusion gene STAR index genome directory', required=True)
