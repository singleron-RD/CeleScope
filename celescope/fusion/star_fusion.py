import celescope.tools.utils as utils
from celescope.tools.star_mixin import Star_mixin, get_opts_star_mixin


class Star_fusion(Star_mixin):
    def __init__(self, args, display_title=None):
        args.genomeDir = args.fusion_genomeDir
        super().__init__(args, display_title)


@utils.add_log
def star_fusion(args):
    with Star_fusion(args) as runner:
        runner.run()


def get_opts_star_fusion(parser, sub_program):
    get_opts_star_mixin(parser, sub_program)
    # will cause `conflicting option string: --genomeDir`
    # parser.add_argument('--genomeDir', help=argparse.SUPPRESS)
    parser.add_argument('--fusion_genomeDir', help='fusion gene STAR index genome directory', required=True)
