from celescope.tools import utils
from celescope.__init__ import __VERSION__
from celescope.chemistry_dict import chemistry_dict


from celescope.tools.step import Step, s_common

from celescope.tools.parse_chemistry import get_chemistry, invalid_debug


class Sample(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.fq1_list = args.fq1.split(",")
        # out
        self.invalid_debug_file = f"{self.out_prefix}_invalid_debug.html"

    @utils.add_log
    def run(self):
        chemistry = get_chemistry(self.assay, self.args.chemistry, self.fq1_list)

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
            help_info='For more information, see <a href="https://github.com/singleron-RD/CeleScope/blob/master/doc/chemistry.md">here</a>',
        )
        self.add_metric(
            name="Software Version",
            value=__VERSION__,
        )
        if self.args.debug:
            invalid_debug(
                chemistry,
                self.fq1_list,
                self.invalid_debug_file,
                pattern=self.args.pattern,
                whitelist=self.args.whitelist,
                linker=self.args.linker,
                use_read=self.args.use_read,
                skip_read=self.args.skip_read,
            )


@utils.add_log
def sample(args):
    with Sample(args) as runner:
        runner.run()


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--fq1", help="read1 fq file")
        parser.add_argument("--linker")
        parser.add_argument(
            "--use_read", help="number of reads used for debug", default=1000, type=int
        )
        parser.add_argument(
            "--skip_read",
            help="Number of reads to skip for debugging. The R1 FASTQ contains N bases in the first hundred reads, so it is better to skip them during debugging",
            default=10000,
            type=int,
        )
    parser.add_argument(
        "--chemistry",
        choices=list(chemistry_dict.keys()),
        help="chemistry version",
        default="auto",
    )
    parser.add_argument("--pattern")
    parser.add_argument(
        "--whitelist",
    )
    return parser
