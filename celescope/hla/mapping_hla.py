import os
from celescope.tools.utils import format_number, log, read_barcode_file


def sub_typing(sub_fq, res_outdir):
    pass


@log
def split_fq(fq, fq_outdir, barcodes):
    pass


@log
def hla_typing(fq_outdir, barcodes, thread):
    pass


@log
def mapping_hla(args):

    sample = args.sample
    outdir = args.outdir
    fq = args.fq
    thread = args.thread
    match_dir = args.match_dir

    # process args
    barcodes, nCell = read_barcode_file(match_dir)

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # split fq
    fq_outdir = f'{outdir}/{split_fq}'
    split_fq(fq, fq_outdir, barcodes)

    # typing
    hla_typing(fq_outdir, barcodes, thread)


def get_opts_mapping_hla(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument("--match_dir", help="match scRNA-Seq dir", required=True)
    parser.add_argument("--thread", help='number of thread', default=1)