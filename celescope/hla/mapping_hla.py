import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
import pysam

import celescope
from celescope.tools.utils import add_log, get_barcode_from_match_dir


@add_log
def razer(fq, outdir, sample, thread):

    # get ref
    root_path = os.path.dirname(celescope.__file__)
    ref = f'{root_path}/data/HLA/hla_reference_rna.fasta'

    # mkdir
    out_bam_dir = f'{outdir}/bam/'
    if not os.path.exists(out_bam_dir):
        os.mkdir(out_bam_dir)

    # out_bam
    out_bam = f'{outdir}/bam/{sample}.bam'

    # run
    cmd = (
        f'razers3 -i 97 -m 99999 --distance-range 0 -pa '
        f'-tc {thread} '
        f'-o {out_bam} '
        f'{ref} '
        f'{fq} '
    )
    razer.logger.info(cmd)
    os.system(cmd)
    return out_bam


@add_log
def split_bam(out_bam, barcodes, outdir, sample):
    '''
    input:
        out_bam: from razers3
        barcodes: cell barcodes
    ouput:
        bam_dict: assign reads to cell barcodes and UMI
        count_dict: UMI counts per cell
        index: assign index(1-based) to cells
    '''

    # init
    count_dict = defaultdict(dict)
    bam_dict = defaultdict(dict)
    index_dict = defaultdict(dict)
    cells_dir = f'{outdir}/cells/'

    # read bam and split
    split_bam.logger.info('reading bam...')
    samfile = pysam.AlignmentFile(out_bam, "rb")
    header = samfile.header
    for read in samfile:
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if barcode in barcodes:
            # keep one read for each UMI
            if umi not in bam_dict[barcode]:
                bam_dict[barcode][umi] = read
            if umi in count_dict[barcode]:
                count_dict[barcode][umi] += 1
            else:
                count_dict[barcode][umi] = 1

    split_bam.logger.info('writing cell bam...')
    # write new bam
    index = 0
    for barcode in barcodes:
        # init
        index += 1
        index_dict[index]['barcode'] = barcode
        index_dict[index]['valid'] = False

        # out
        if barcode in bam_dict:
            cell_dir = f'{cells_dir}/cell{index}'
            cell_bam_file = f'{cell_dir}/cell{index}.bam'
            if not os.path.exists(cell_dir):
                os.makedirs(cell_dir)
            index_dict[index]['valid'] = True
            cell_bam = pysam.AlignmentFile(
                f'{cell_bam_file}', "wb", header=header)
            for umi in bam_dict[barcode]:
                read = bam_dict[barcode][umi]
                cell_bam.write(read)
            cell_bam.close()

    # out df_index
    df_index = pd.DataFrame(index_dict).T
    df_index.index.name = 'cell_index'
    index_file = f'{outdir}/{sample}_cell_index.tsv'
    df_index.to_csv(index_file, sep='\t')

    # out count_dict
    df_count = pd.DataFrame(count_dict).T
    df_temp = df_count.stack()
    df_temp = df_temp.astype(int)
    df_temp = df_temp.reset_index()
    df_temp.columns = ['barcode', 'UMI', 'read_count']
    count_file = f'{outdir}/{sample}_UMI_count.tsv'
    df_temp.to_csv(count_file, sep='\t', index=False)

    return index_file, count_file


def sub_typing(bam):

    outdir = os.path.dirname(bam)
    prefix = os.path.basename(bam).strip('.bam')
    cmd = (
        f'OptiTypePipeline.py '
        f'--input {bam} '
        f'--rna '
        f'--outdir {outdir} '
        f'--prefix {prefix} '
        f'>/dev/null 2>&1 '
    )
    os.system(cmd)


def read_index(index_file):
    df_index = pd.read_csv(index_file, sep='\t', index_col=0, dtype=object)
    df_valid = df_index[df_index['valid'] == 'True']
    return df_valid


@add_log
def hla_typing(index_file, outdir, thread):
    all_res = []
    df_valid = read_index(index_file)

    bams = [
        f'{outdir}/cells/cell{index}/cell{index}.bam' for index in df_valid.index]
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(sub_typing, bams):
            all_res.append(res)


@add_log
def summary(index_file, outdir, sample):

    n = 0
    df_valid = read_index(index_file)

    for index in df_valid.index:
        try:
            sub_df = pd.read_csv(
                f'{outdir}/cells/cell{index}/cell{index}_result.tsv', sep='\t', index_col=0)
        except FileNotFoundError:
            continue
        n += 1
        sub_df['barcode'] = df_valid.loc[index, :]['barcode']
        sub_df['cell_index'] = index
        if n == 1:
            all_df = sub_df
        else:
            all_df = all_df.append(sub_df, ignore_index=True)
    all_df['Reads'] = all_df['Reads'].apply(int)
    all_df = all_df[all_df['Reads'] != 0]
    all_df = all_df.drop('Objective', axis=1)
    out_file = f'{outdir}/{sample}_typing.tsv'
    all_df.to_csv(out_file, sep='\t', index=False)


@add_log
def mapping_hla(args):

    sample = args.sample
    outdir = args.outdir
    fq = args.fq
    thread = int(args.thread)
    match_dir = args.match_dir

    # process args
    barcodes, _nCell = get_barcode_from_match_dir(match_dir)

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # razer
    out_bam = razer(fq, outdir, sample, thread)

    # split bam
    index_file, _count_file = split_bam(out_bam, barcodes, outdir, sample)

    # typing
    hla_typing(index_file, outdir, thread)

    # summary
    summary(index_file, outdir, sample)


def get_opts_mapping_hla(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        "--match_dir", help="match scRNA-Seq dir", required=True)
    parser.add_argument("--thread", help='number of thread', default=1)
