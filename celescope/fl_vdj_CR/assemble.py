import os
from celescope.fl_vdj_CR.VDJ_Mixin import VDJ_Mixin, get_opts_VDJ_Mixin
from celescope.tools.step import s_common


class Assemble(VDJ_Mixin):
    """
    ## Features

    - TCR/BCR Assemble.

    ## Output
    - `03.assemble/{sample}/outs/` Recording assemble results.

    """
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.fqs_dir = os.path.abspath(args.fqs_dir)

    def run(self):
        self.run_assemble()


def assemble(args):
    assemble_obj = Assemble(args)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    get_opts_VDJ_Mixin(parser)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir after convert', required=True)
    return parser

