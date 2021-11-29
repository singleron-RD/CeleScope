from celescope.tools.utils import STAR_util, add_log


@add_log
def mapping_mut(args):
    STAR_util(
        args.sample,
        args.outdir,
        args.fq,
        args.indel_genomeDir,
        args.thread,
        args.outFilterMatchNmin,
        outBAMsortingBinsN=200,
    )


def get_opts_mapping_mut(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--indel_genomeDir',
        help='insertion or deletion STAR indexed genome directory',
        required=True)
    parser.add_argument("--thread", help='STAR thread', default=1)
    parser.add_argument("--outFilterMatchNmin", help='STAR outFilterMatchNmin', default=35)
