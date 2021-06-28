import os

import pandas as pd

import celescope.tools.utils as utils
from celescope.__init__ import __VERSION__, ASSAY_DICT
from celescope.tools.__init__ import __PATTERN_DICT__
from celescope.tools.barcode import Chemistry
from celescope.tools.step import Step, s_common


@utils.add_log
def sample(args):

    step_name = "sample"
    step = Step(args, step_name)

    sample_name = args.sample
    assay = args.assay
    assay_description = ASSAY_DICT[assay]
    version = __VERSION__
    outdir = args.outdir
    chemistry = args.chemistry

    # get chemistry
    if chemistry == 'auto':
        fq1 = args.fq1
        ch = Chemistry(fq1)
        chemistry = ch.check_chemistry()
        chemistry = ",".join(set(chemistry))
    else:
        chemistry = args.chemistry

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    stat = pd.DataFrame({
        "item": ["Sample ID", "Assay", "Chemistry", "Software Version"],
        "count": [sample_name, assay_description, chemistry, version],
    },
        columns=["item", "count"]
    )
    stat_file = outdir + "/stat.txt"
    stat.to_csv(stat_file, sep=":", header=None, index=False)

    step.clean_up()

    return chemistry


def get_opts_sample(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='read1 fq file')
    parser.add_argument('--chemistry', choices=list(__PATTERN_DICT__.keys()), help='chemistry version', default='auto')
    return parser
