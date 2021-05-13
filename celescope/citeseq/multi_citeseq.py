import os
from collections import defaultdict
from .__init__ import __STEPS__, __ASSAY__
from celescope.tools.utils import *


def main():

    # init
    assay = __ASSAY__
    steps = __STEPS__
    conda = os.environ['CONDA_DEFAULT_ENV']
    app = 'celescope'

    # parser
    parser = multi_opts(assay)
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument("--fq_pattern", help="tag read2 pattern", required=True)
    parser.add_argument("--linker_fasta", help="linker fasta")
    parser.add_argument("--barcode_fasta", help="barcode fasta", required=True)
    args = parser.parse_args()

    # read args
    outdir = args.outdir
    chemistry = args.chemistry
    pattern = args.pattern
    whitelist = args.whitelist
    linker = args.linker
    lowQual = args.lowQual
    lowNum = args.lowNum
    mod = args.mod
    rm_files = args.rm_files
    minimum_length = args.minimum_length

    # parse mapfile
    fq_dict, match_dict = parse_map_col4(args.mapfile, None)

    # link
    link_data(outdir, fq_dict)

    # custom args
    thread = args.thread
    fq_pattern = args.fq_pattern
    linker_fasta = args.linker_fasta
    barcode_fasta = args.barcode_fasta

    # mk log dir
    logdir = outdir + '/log'
    os.system('mkdir -p %s' % (logdir))

    # script init
    sjm_cmd = 'log_dir %s\n' % (logdir)
    sjm_order = ''
    shell_dict = defaultdict(str)

