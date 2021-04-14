import os
from celescope.tools.utils import *
from celescope.tools.STAR import Step_mapping


@add_log
def STAR_fusion(args):
    s = Step_mapping(
        args.sample,
        args.outdir,
        args.assay,
        args.thread,
        args.input_read,
        args.genomeDir,
        outFilterMatchNmin = args.outFilterMatchNmin,
        sort_BAM=False
    )
    s.STAR()
    s.sort_bam()


def get_opts_STAR_fusion(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument("--input_read", required=True)

    parser.add_argument(
        '--genomeDir',
        help='fusion gene STAR index genome directory',
        required=True)
    parser.add_argument("--outFilterMatchNmin", help='STAR outFilterMatchNmin', default=35)
    parser.add_argument('--starMem', help='starMem', default=30)
