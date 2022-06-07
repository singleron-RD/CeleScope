from celescope.fl_vdj_CR.VDJ_Mixin import VDJ_Mixin, get_opts_VDJ_Mixin
from celescope.tools.step import s_common


class Convert(VDJ_Mixin):
    """
    Features

    - Format barcodes and UMIs.

    Output        
    - `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

    - `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads.

    - `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads.
    """
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        self.fq2 = args.fq2

    def run(self):
        self.run_convert()


def convert(args):
    step_name = 'convert'
    convert_obj = Convert(args, step_name)
    convert_obj.run()


def get_opts_convert(parser, sub_program):
    get_opts_VDJ_Mixin(parser)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fq2', help='R2 read file', required=True)
    return parser
