import argparse
from celescope.tools.utils import *
from .Mapping_tag import Mapping_tag

def get_opts_mapping_tag(parser, sub_program):
    parser.add_argument("--fq_pattern", help="read2 fastq pattern")
    parser.add_argument("--linker_fasta", help="linker fasta")
    parser.add_argument("--barcode_fasta", help="barcode fasta")
    if sub_program:
        s_common(parser)
        parser.add_argument("--fq", help="clean read2", required=True)


def mapping_tag(args):
    sample = args.sample
    outdir = args.outdir
    assay = args.assay
    fq = args.fq
    fq_pattern = args.fq_pattern
    linker_fasta = args.linker_fasta
    barcode_fasta = args.barcode_fasta

    mapping_tag_object = Mapping_tag(
        sample,
        outdir,
        assay,
        fq,
        fq_pattern,
        linker_fasta,
        barcode_fasta,
    )
    mapping_tag_object.run()
