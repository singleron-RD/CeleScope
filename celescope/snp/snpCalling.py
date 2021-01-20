import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import celescope
import pysam
import numpy as np
import pandas as pd
import logging
from celescope.tools.utils import format_number, log, read_barcode_file
from celescope.tools.utils import format_stat
from celescope.tools.utils import read_one_col, gene_convert, glob_genomeDir
from celescope.tools.report import reporter
from mutract.utils import Mutract


@log
def snpCalling(args):

    sample = args.sample
    outdir = args.outdir
    thread = int(args.thread)
    match_dir = args.match_dir
    bam = args.bam
    genomeDir = args.genomeDir
    gene_list_file = args.gene_list

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # get args
    _refFlat, _gtf, fasta = glob_genomeDir(genomeDir, fa=True)
    _match_barcode, _cell_total, match_barcode_file = read_barcode_file(match_dir, return_file=True)

    # mutract
    obj = Mutract(
        outdir, sample, bam, fasta, 
        match_barcode_file, thread=thread, gene_file=gene_list_file
    )
    obj.run()


def get_opts_snpCalling(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--assay', help='assay', required=True)
        parser.add_argument("--thread", help='number of thread', default=1)
        parser.add_argument("--bam", help='featureCounts bam', required=True)
        parser.add_argument("--genomeDir", help='genomeDir', required=True)
    parser.add_argument(
        "--match_dir", help="match scRNA-Seq dir", required=True)
    parser.add_argument("--gene_list", help='gene_list', required=True)
