import os
import subprocess
from celescope.tools.utils import log
from celescope.tools.STAR import Step_mapping

@log
def STAR_virus(args):

    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % (args.outdir))

    sample = f'{args.sample}_virus'
    mapping = Step_mapping(
        sample, 
        args.outdir, 
        args.assay, 
        args.thread,
        args.input_read, 
        args.virus_genomeDir,
        outFilterMultimapNmax=10,
        outFilterMatchNmin=args.outFilterMatchNmin,
        sort_BAM=False
    )
    mapping.STAR()
    mapping.sort_bam()
    mapping.index_bam()


def get_opts_STAR_virus(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--input_read", required=True)
        parser.add_argument("--thread", help='STAR thread', default=1)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--virus_genomeDir',
        help='virus genome dir',
        required=True)
    parser.add_argument("--outFilterMatchNmin", help='STAR outFilterMatchNmin', default=35)
