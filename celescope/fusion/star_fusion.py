from celescope.tools import utils
from celescope.tools.star_mixin import Star_mixin, get_opts_star_mixin


class Star_fusion(Star_mixin):
    """
    ## Features
    - The reads were aligned to the known fusion sites(specified by `--fusion_genomeDir`) using STAR.
    Please note that [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) is not used here.
    """

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
    parser.add_argument(
        "--fusion_genomeDir",
        help="fusion gene STAR index genome directory",
        required=True,
    )
