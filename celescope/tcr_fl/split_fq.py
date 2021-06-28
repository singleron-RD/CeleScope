import os
from collections import defaultdict

import pandas as pd
import pysam

from celescope.tcr_fl.barcode_index import Barcode_index
from celescope.tools.utils import add_log, genDict, read_barcode_file


@add_log
def get_nCell_barcodes(fq, nCell):
    '''
    get top nCell's barcodes(rank by UMI counts)
    '''
    count_dict = genDict(dim=2)
    barcode_dict = {}
    with pysam.FastxFile(fq) as fh:
        for entry in fh:
            attr = entry.name.split('_')
            barcode = attr[0]
            umi = attr[1]
            count_dict[barcode][umi] += 1
    for barcode in count_dict:
        barcode_dict[barcode] = len(count_dict[barcode])
    barcodes = pd.DataFrame.from_dict(barcode_dict, orient='index').sort_values(
        0, ascending=False).iloc[0:nCell, ].index
    return barcodes


@add_log
def split_run(fq, fq_outdir, barcodes=None, nCell=None):
    '''
    split fastq 
    '''
    if not os.path.exists(fq_outdir):
        os.makedirs(fq_outdir)
    if nCell and nCell != 'None':
        barcodes = get_nCell_barcodes(fq, nCell)
    bi = Barcode_index(barcodes)
    entry_dict = defaultdict(list)
    with pysam.FastxFile(fq) as fh:
        for entry in fh:
            attr = entry.name.split('_')
            barcode = attr[0]
            if barcode in barcodes:
                cell_index = bi.index_dict[barcode]
                entry_dict[cell_index].append(entry)

    # write to file
    for cell_index in entry_dict:
        with open(f'{fq_outdir}/{cell_index}.fq', 'w') as f:
            for entry in entry_dict[cell_index]:
                f.write(str(entry) + '\n')
    return bi


@add_log
def split_fq(args):
    nCell = args.nCell
    outdir = args.outdir
    sample = args.sample
    match_dir = args.match_dir

    if match_dir and match_dir != 'None':
        barcodes, _nCell = read_barcode_file(args.match_dir)
    else:
        barcodes = ''
    fq_outdir = f'{args.outdir}/fastq'
    if nCell and nCell != 'None':
        nCell = int(nCell)
    bi = split_run(args.fq, fq_outdir, barcodes, nCell)
    index_file = f'{outdir}/{sample}_index.tsv'
    bi.df_index.to_csv(index_file, sep='\t')


def get_opts_split_fq(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        "--match_dir", help="match scRNA-Seq dir")
    parser.add_argument("--nCell", help="select top N cell")
