from celescope.tools import utils
from celescope.__init__ import __VERSION__
from celescope.tools.__init__ import PATTERN_DICT

from celescope.tools.step import Step, s_common

from sccore.parse_protocol import AutoRNA


def parse_chemistry(chemistry, fq1_list):
    if chemistry == "auto":
        return AutoRNA(fq1_list).get_protocol()
    return chemistry


class Sample(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.fq1_list = args.fq1.split(",")

    @utils.add_log
    def run(self):
        chemistry = parse_chemistry(self.args.chemistry, self.fq1_list)

        if chemistry == "bulk_rna":
            chemistry = f"accuracode{self.args.wells}"

        self.add_metric(
            name="Sample ID",
            value=self.sample,
        )
        self.add_metric(
            name="Assay",
            value=self.assay,
        )
        self.add_metric(
            name="Chemistry",
            value=chemistry,
            help_info='For more information, see <a href="https://github.com/singleron-RD/CeleScope/blob/cbf5df39b74628fbcc1d9265728dcfe86b6989f2/doc/chemistry.md">here</a>',
        )
        self.add_metric(
            name="Software Version",
            value=__VERSION__,
        )


@utils.add_log
def sample(args):
    with Sample(args) as runner:
        runner.run()


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--fq1", help="read1 fq file")
    parser.add_argument(
        "--chemistry",
        choices=list(PATTERN_DICT.keys()),
        help="chemistry version",
        default="auto",
    )
    parser.add_argument(
        "--wells",
        type=int,
        help="The AccuraCode wells used (384 or 96), default 384.",
        default=384,
    )
    return parser
