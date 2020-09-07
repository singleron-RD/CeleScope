import os
import gzip
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import celescope
import pysam
import pandas as pd
from celescope.tools.utils import format_number, log, read_barcode_file


@log
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


@log
def split_bam(out_bam, outdir, barcodes, sample):

    # init
    count_dict = defaultdict(dict)
    bam_dict = defaultdict(list)
    index_dict = defaultdict(dict)
    cells_dir = f'{outdir}/cells/'

    # read bam and split
    samfile = pysam.AlignmentFile(out_bam, "rb")
    header = samfile.header
    for read in samfile:
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if barcode in barcodes:
            bam_dict[barcode].append(read)
            if umi in count_dict[barcode]:
                count_dict[barcode][umi] += 1
            else:
                count_dict[barcode][umi] = 1
    
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
            for read in bam_dict[barcode]:
                cell_bam.write(read)
            cell_bam.close()

    # out df_index
    df_index = pd.DataFrame(index_dict).T
    index_file = f'{outdir}/{sample}_cell_index.tsv'
    df_index.to_csv(index_file, sep='\t')

    # out count_dict
    df_count = pd.DataFrame(count_dict).T
    count_file = f'{outdir}/{sample}_UMI_count.tsv'
    df_count.to_csv(count_file, sep='\t')

    return count_dict, index_dict


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


@log
def hla_typing(index_dict, outdir, thread):
    all_res = []
    valid_index_dict = {}
    for index in index_dict:
        if index_dict[index]['valid']:
            valid_index_dict[index] = index_dict[index]['barcode']

    bams = [f'{outdir}/cells/cell{index}/cell{index}.bam' for index in valid_index_dict]
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(sub_typing, bams):
            all_res.append(res)
    return valid_index_dict


@log
def summary(valid_index_dict, outdir):
    pass


@log
def mapping_hla(args):

    sample = args.sample
    outdir = args.outdir
    fq = args.fq
    thread = int(args.thread)
    match_dir = args.match_dir

    # process args
    barcodes, _nCell = read_barcode_file(match_dir)

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # razer
    out_bam = razer(fq, outdir, sample, thread)

    # split bam
    _count_dict, index_dict = split_bam(out_bam, outdir, barcodes, sample)

    # typing
    valid_index_dict = hla_typing(index_dict, outdir, thread)

    # summary
    summary()



def get_opts_mapping_hla(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument("--match_dir", help="match scRNA-Seq dir", required=True)
    parser.add_argument("--thread", help='number of thread', default=1)