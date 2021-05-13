import os
from collections import defaultdict
from .__init__ import __ASSAY__
from celescope.tools.utils import link_data, multi_opts, parse_map_col4


def main():

    # init
    assay = __ASSAY__
    os.environ['CONDA_DEFAULT_ENV']

    # parser
    parser = multi_opts(assay)
    parser.add_argument('--thread', help='thread', default=6)
    parser.add_argument("--fq_pattern", help="tag read2 pattern", required=True)
    parser.add_argument("--linker_fasta", help="linker fasta")
    parser.add_argument("--barcode_fasta", help="barcode fasta", required=True)
    args = parser.parse_args()

    # read args
    outdir = args.outdir
    args.chemistry
    args.pattern
    args.whitelist
    args.linker
    args.lowQual
    args.lowNum
    args.mod
    args.rm_files
    args.minimum_length

    # parse mapfile
    fq_dict, match_dict = parse_map_col4(args.mapfile, None)

    # link
    link_data(outdir, fq_dict)

    # custom args
    args.thread
    args.fq_pattern
    args.linker_fasta
    args.barcode_fasta

    # mk log dir
    logdir = outdir + '/log'
    os.system('mkdir -p %s' % (logdir))

    # script init
    sjm_cmd = 'log_dir %s\n' % (logdir)
    defaultdict(str)

