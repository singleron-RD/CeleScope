

import celescope.tools.utils as utils
from celescope.__init__ import __VERSION__
from celescope.tools.__init__ import __PATTERN_DICT__
from celescope.tools.barcode import Chemistry
from celescope.tools.step import Step, s_common


class Sample(Step):
    def __init__(self, args):
        Step.__init__(self, args)
        self.assay_description = utils.get_assay_text(self.assay)
        self.version = __VERSION__
        self.chemistry = args.chemistry

    @utils.add_log
    def run(self):
        if self.chemistry == 'auto':
            fq1 = self.args.fq1
            ch = Chemistry(fq1)
            chemistry = ch.check_chemistry()
            chemistry = ",".join(set(chemistry))
        else:
            chemistry = self.chemistry

        self.add_metric(
            name='Sample ID',
            value=self.sample,
        )
        self.add_metric(
            name='Assay',
            value=self.assay_description,
        )
        self.add_metric(
            name='Chemistry',
            value=chemistry,
            help_info='Chemistry of the input fastqs',
        )
        self.add_metric(
            name='Software Version',
            value=self.version,
        )


@utils.add_log
def sample(args):

    with Sample(args) as runner:
        runner.run()


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='read1 fq file')
    parser.add_argument('--chemistry', choices=list(__PATTERN_DICT__.keys()), help='chemistry version', default='auto')
    return parser
