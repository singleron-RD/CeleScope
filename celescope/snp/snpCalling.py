import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import celescope
import pysam
import pandas as pd
from celescope.tools.utils import format_number, log, read_barcode_file
from celescope.tools.utils import read_one_col, gene_convert, glob_genomeDir


@log
def split_bam(bam, barcodes, outdir, sample, gene_id_name_dic, min_query_length):
    '''
    input:
        bam: bam from feauturCounts
        barcodes: cell barcodes
        gene_id_name_dic: id name dic
        min_query_length: minimum query length

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
    samfile = pysam.AlignmentFile(bam, "rb")
    header = samfile.header
    for read in samfile:
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        if not read.has_tag('XT'):
            continue
        gene = read.get_tag('XT')
        query_length = read.infer_query_length()
        if (barcode in barcodes) and (gene in gene_id_name_dic) and (query_length >= min_query_length):
            gene_name = gene_id_name_dic[gene]
            # keep one read for each UMI
            if umi not in bam_dict[barcode]:
                bam_dict[barcode][umi] = read
            # count
            if gene_name not in count_dict[barcode]:
                count_dict[barcode][gene_name] = {}
            if umi in count_dict[barcode][gene_name]:
                count_dict[barcode][gene_name][umi] += 1
            else:
                count_dict[barcode][gene_name][umi] = 1

    split_bam.logger.info('writing cell bam...')
    # write new bam
    index = 0
    for barcode in barcodes:
        # init
        index += 1
        index_dict[index]['barcode'] = barcode
        index_dict[index]['valid'] = False

        # out bam
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
    df_count = pd.DataFrame(columns=['barcode', 'gene', 'UMI', 'read_count'])
    for barcode in count_dict:
        for gene in count_dict[barcode]:
            for umi in count_dict[barcode][gene]:
                read_count = count_dict[barcode][gene][umi]
                df_count = df_count.append({
                    'barcode': barcode,
                    'gene': gene,
                    'UMI': umi,
                    'read_count': read_count,
                }, ignore_index=True)
    count_file = f'{outdir}/{sample}_count.tsv'
    df_count.to_csv(count_file, sep='\t', index=False)

    return index_file, count_file


def call_snp(index, outdir, fasta):

    bam = f'{outdir}/cells/cell{index}/cell{index}.bam'
    # sort
    sorted_bam = f'{outdir}/cells/cell{index}/cell{index}_sorted.bam'
    cmd_sort = (
        f'samtools sort {bam} -o {sorted_bam}'
    )
    os.system(cmd_sort)

    # mpileup
    out_bcf = f'{outdir}/cells/cell{index}/cell{index}.bcf'
    cmd_mpileup = (
        f'samtools mpileup -g -f '
        f'{fasta} '
        f'{sorted_bam} | '
        f'bcftools call -mv -Ob '
        f'-o {out_bcf}'
    )
    os.system(cmd_mpileup)

    # view
    out_vcf = f'{outdir}/cells/cell{index}/cell{index}.vcf'
    cmd_view = (
        f'bcftools view {out_bcf} > {out_vcf}'
    )
    os.system(cmd_view)


def read_index(index_file):
    df_index = pd.read_csv(index_file, sep='\t', index_col=0, dtype=object)
    df_valid = df_index[df_index['valid'] == 'True']
    return df_valid


@log
def call_all_snp(index_file, outdir, thread, fasta):
    all_res = []
    df_valid = read_index(index_file)
    index_arg = df_valid.index
    outdir_arg = [outdir] * len(index_arg)
    fasta_arg = [fasta] * len(index_arg)
    with ProcessPoolExecutor(thread) as pool:
        for res in pool.map(call_snp, index_arg, outdir_arg, fasta_arg):
            all_res.append(res)


@log
def summary(index_file, outdir, sample):
    
    n = 0
    df_valid = read_index(index_file)
    
    for index in df_valid.index:
        try:
            sub_df = pd.read_csv(
                f'{outdir}/cells/cell{index}/cell{index}_result.tsv', sep='\t', index_col=0)
        except Exception:
            continue
        n += 1
        sub_df['barcode'] = df_valid.loc[index, :]['barcode']
        sub_df['cell_index'] = index
        if n == 1:
            all_df = sub_df
        else:
            all_df = all_df.append(sub_df, ignore_index=True)
    all_df['Reads'] = all_df['Reads'].apply(lambda x: int(x))
    all_df = all_df[all_df['Reads'] != 0]
    all_df = all_df.drop('Objective', axis=1)
    out_file = f'{outdir}/{sample}_typing.tsv'
    all_df.to_csv(out_file, sep='\t', index=False)


@log
def convert(gene_list_file, gtf):
    gene_list_name, _count = read_one_col(gene_list_file)
    id_name = gene_convert(gtf)
    name_id = {}
    for id in id_name:
        name = id_name[id]
        name_id[name] = id
    gene_id_name_dic = {}
    for gene_name in gene_list_name:
        gene_id = name_id[gene_name]
        gene_id_name_dic[gene_id] = gene_name
    return gene_id_name_dic


@log
def snpCalling(args):

    sample = args.sample
    outdir = args.outdir
    thread = int(args.thread)
    match_dir = args.match_dir
    bam = args.bam
    genomeDir = args.genomeDir
    gene_list_file = args.gene_list
    min_query_length = args.min_query_length

    # process args
    barcodes, _nCell = read_barcode_file(match_dir)

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    # get genome file
    _refFlat, gtf, fasta = glob_genomeDir(genomeDir, fa=True)

    # convert gene
    gene_id_name_dic = convert(gene_list_file, gtf)

    # split bam
    index_file, _count_file = split_bam(
        bam, barcodes, outdir, sample, gene_id_name_dic, min_query_length)

    # snp
    call_all_snp(index_file, outdir, thread, fasta)

    # summary
    # summary(index_file, outdir, sample)


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
    parser.add_argument("--min_query_length", help='minimum query length', default=35)
