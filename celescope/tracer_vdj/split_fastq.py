import pysam
from collections import defaultdict
import os
import argparse
import datetime
import pandas as pd
from Bio.Seq import Seq
import glob
from celescope.tools import utils
from celescope.tools.utils import *


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
    with open(f'{barcodes_path}/reversed_barcodes.tsv', 'w') as fh:
        for barcode in barcodes:
            barcode = Seq(barcode)
            barcode_reversed = barcode.reverse_complement()
            bc = str(barcode_reversed)
            fh.write(bc + '\n')

    with open(f'{barcodes_path}/reversed_barcodes.tsv') as res:
        res = res.readlines()
    return res


@utils.add_log
def get_fastq_to_assemble(fq_outdir, fq, barcodes):
    """
    split_fastq
    """
    if not os.path.exists(fq_outdir):
        os.makedirs(fq_outdir)
    
    barcode_reads_dict = defaultdict(list)  # all barcodes from BCR vdj_dir paired with reads
    umi_count = defaultdict()
    reads_count_dict = {}  # all barcodes and reads num for each barcode
    all_barcodes = []  # all barcodes
    with pysam.FastxFile(fq) as fq:
        for entry in fq:
            attr = entry.name.split('_')
            barcode = attr[0]
            umi = attr[1]
            umi_count[barcode][umi] += 1
            all_barcodes.append(barcode)
            barcode_reads_dict[barcode].append(entry)
        for barcode in list(barcode_reads_dict.keys()):
            reads_count_dict[barcode] = len(barcode_reads_dict[barcode])
    
        umi_count_df = pd.DataFrame([(k, list(v.keys())[0], list(v.values())[0]) for k, v in umi_count.items()], columns=['Barcode', 'umi', 'umi_reads_count'])

        umi_df = umi_count_df.groupby(['Barcode']).agg({'UMI': 'count'})

        umi_df.to_csv(f'{fq_outdir}/../umi_count.tsv', sep='\t')
        
        barcodes_for_match = []
        for barcode in barcodes:
                barcode = barcode.strip('\n')
                barcodes_for_match.append(barcode)
        barcodes_to_use = list(set(barcodes_for_match).intersection(set(all_barcodes)))
            # barcodes in both RNA data and BCR data

    barcode_reads_useful = {barcode: barcode_reads_dict[barcode] for barcode in barcodes_to_use}


    barcodes_reads_count = {barcode: reads_count_dict[barcode] for barcode in
                            list(barcode_reads_useful.keys())}

    barcodes_reads_cal = pd.DataFrame.from_dict(barcodes_reads_count, orient='index',columns=['counts'])
    barcodes_reads_cal = barcodes_reads_cal.reset_index().rename(columns={'index': 'barcode'})
    barcodes_reads_cal = barcodes_reads_cal.sort_values(by='counts', ascending=False)

    barcodes_reads_cal.to_csv(f'{fq_outdir}/../reads_count.tsv', sep='\t')

    stat_string = 'All_cells:\t{}\nmatched_cell:\t{}\n'.format(len(all_barcodes), len(barcode_reads_useful))
    with open(f'{fq_outdir}/../stat.txt', 'w') as s:
        s.write(stat_string)

    i = 1
    for barcode in list(barcode_reads_useful.keys()):

        with open(f'{fq_outdir}/{i}.fq', 'w') as f:
            for entry in barcode_reads_useful[barcode]:
                f.write(str(entry) + '\n')
        if i % 1000 == 0:
            get_fastq_to_assemble.logger.info(f'processed {i} cells')

        if i == len(list(barcode_reads_useful.keys())):
            get_fastq_to_assemble.logger.info(f'finnaly get {i} cells')

        i += 1
        


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
    


