

from celescope.tools import utils
from celescope.__init__ import __VERSION__
from celescope.tools.__init__ import PATTERN_DICT
from celescope.tools.barcode import Chemistry
from celescope.tools.step import Step, s_common


def add_kit_version(chemistry):
    kit_dict = {
        '1': 'no longer in use',
        '2': 'kit V1',
        '3': 'kit V2',
    }
    if chemistry.startswith('scopeV'):
        s = chemistry.replace('scopeV', '')
        chem_version = s[0]
        kit = kit_dict[chem_version]
        chemistry = f'{chemistry} ({kit})'
    
    return chemistry

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
            display=add_kit_version(chemistry),
            help_info='For more information, see <a href="https://github.com/singleron-RD/CeleScope/blob/master/docs/chemistry.md">here</a>',
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
    parser.add_argument('--chemistry', choices=list(PATTERN_DICT.keys()), help='chemistry version', default='auto')
    return parser
