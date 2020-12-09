import os
from celescope.tools.utils import log
from celescope.tools.STAR import Step_mapping


@log
def STAR_fusion(args):
    s = Step_mapping(
        args.sample,
        args.outdir,
        args.assay,
        args.thread,
        args.input_read,
        args.genomeDir,
        outFilterMatchNmin = args.outFilterMatchNmin,
    )
    s.STAR()


def get_opts_STAR_fusion(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--input_read", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--genomeDir',
        help='fusion gene STAR index genome directory',
        required=True)
    parser.add_argument("--thread", help='STAR thread', default=1)
    parser.add_argument("--outFilterMatchNmin", help='STAR outFilterMatchNmin', default=35)
