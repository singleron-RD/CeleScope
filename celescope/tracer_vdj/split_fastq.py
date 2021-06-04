import pysam
from collections import defaultdict
import os
import argparse
import datetime
import pandas as pd
from Bio.Seq import Seq
import glob
from celescope.tools import utils
from celescope.tools.Step import Step, s_common


@utils.add_log
def annotation_barcodes(match_dir, type):
    
    cluster_data = glob.glob(f'{match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')
    cluster_data = cluster_data[0]
    cluster_type = pd.read_csv(cluster_data, sep='\t')

    # filter barcodes
    if type == 'TCR':
        clusters = list(cluster_type[cluster_type['cell_type'] == 'T cells']['cluster'])
    elif type == 'BCR':
        clusters = list(cluster_type[cluster_type['cell_type'] == 'B cells']['cluster'])

    tsne = glob.glob(f'{match_dir}/06.analysis/*_tsne_coord.tsv')
    tsne = tsne[0]
    tsne_coord = pd.read_csv(tsne, sep='\t', index_col=0)

    barcodes = []
    for cluster in clusters:
        tmp = tsne_coord[tsne_coord['cluster'] == cluster].index.tolist()
        barcodes += tmp
    # write barcodes
    barcodes_path = glob.glob(f'{match_dir}/06.analysis/*_auto_assign/')
    barcodes_path = barcodes_path[0]

    res = [] 
    with open(f'{barcodes_path}/reversed_barcodes.tsv', 'w') as fh:
        for barcode in barcodes:
            barcode = Seq(barcode)
            barcode_reversed = barcode.reverse_complement()
            bc = str(barcode_reversed)
            res.append(bc)
            fh.write(bc + '\n')

    return res


@utils.add_log
def get_fastq_to_assemble(fq_outdir, fq, barcodes):
    """
    split_fastq
    """
    if not os.path.exists(fq_outdir):
        os.makedirs(fq_outdir)
    
    barcode_reads_dict = defaultdict(list)  # all barcodes from BCR vdj_dir paired with reads
    # umi_count = defaultdict(list)
    reads_count_dict = {}  # all barcodes and reads num for each barcode
    umi_count_dict = defaultdict(list)
    umi_count = {}

    with pysam.FastxFile(fq) as fq:
        for entry in fq:
            attr = entry.name.split('_')
            barcode = attr[0]
            umi = attr[1]
            if barcode in barcodes:
                barcode_reads_dict[barcode].append(entry)
            if umi_count_dict[barcode].count(umi) == 0:
                umi_count_dict[barcode].append(umi)
        for barcode in barcodes:
            reads_count_dict[barcode] = len(barcode_reads_dict[barcode])
      
        for barcode in list(umi_count_dict.keys()):
            umi_count[barcode] = len(umi_count_dict[barcode])

    df_umi = pd.DataFrame.from_dict(umi_count, orient='index',columns=['UMI'])
    df_umi = df_umi.reset_index().rename(columns={'index': 'Barcode'})

    df_umi.to_csv(f'{fq_outdir}/../umi_count.tsv', sep='\t')

    reads_count = pd.DataFrame.from_dict(reads_count_dict, orient='index',columns=['readcount'])
    reads_count = reads_count.reset_index().rename(columns={'index': 'Barcode'})

    df_f = pd.merge(reads_count, df_umi, on='Barcode', how='inner')

    df_f = df_f.set_index('Barcode')

    i = 1

    for barcode in barcodes:

        df_f.loc[barcode, 'cell_name'] = i

        with open(f'{fq_outdir}/{i}.fq', 'w') as f:
            for entry in barcode_reads_dict[barcode]:
                f.write(str(entry) + '\n')

        if i % 1000 == 0:
            get_fastq_to_assemble.logger.info(f'processed {i} cells')

        if i == len(barcodes):
            get_fastq_to_assemble.logger.info(f'finally get {i} cells')

        i += 1
    
    df_f = df_f.astype(int)
    df_f.to_csv(f'{fq_outdir}/../reads_count.tsv', sep='\t')
        


def split_fastq(args):
    type = args.type
    match_dir = args.match_dir
    sample = args.sample
    outdir = args.outdir
    assay = args.assay
    fq = args.fq

    fq_outdir = f'{outdir}/fastq'
    barcodes = annotation_barcodes(match_dir, type)
        
    get_fastq_to_assemble(fq_outdir, fq, barcodes)


def get_opts_split_fastq(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq', required=True)
        parser.add_argument('--match_dir', help='matched rna_dir')
    parser.add_argument('--type', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    


