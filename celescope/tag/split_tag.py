"""
split scRNA-Seq fastq file(01.barcode/{sample}_2.fq)
"""
import glob
import os
import itertools
from collections import defaultdict

import pysam
import pandas as pd

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools.__init__ import FILTERED_MATRIX_DIR_SUFFIX
from celescope.__init__ import HELP_DICT
from celescope.tools.matrix import CountMatrix


def get_clonotypes_table(df):
    chains = sorted(set(df['chain'].tolist()))
    res = pd.DataFrame(columns=['barcode'])
    for c in chains:
        df_c = df[df['chain'] == c][['barcode', 'aaSeqCDR3', 'nSeqCDR3']]
        df_c = df_c.rename(columns={'aaSeqCDR3': f'{c}_aaSeqCDR3', 'nSeqCDR3': f'{c}_nSeqCDR3'})
        res = pd.merge(res, df_c, on='barcode', how='outer').fillna('NaN')
    group_l = res.columns.tolist()
    group_l.remove('barcode')
    clonotypes = res.groupby(group_l, as_index=False).agg({'barcode': 'count'})
    clonotypes = clonotypes.rename(columns={'barcode': 'barcode_count'})
    sum_cb = clonotypes['barcode_count'].sum()
    clonotypes['percent'] = clonotypes['barcode_count'].apply(lambda x: round(x/sum_cb, 2))
    clonotypes['clonetype_ID'] = [i for i in range(1, clonotypes.shape[0]+1)]
    group_l.insert(0, 'clonetype_ID')
    group_l.append('barcode_count')
    group_l.append('percent')
    clonotypes = clonotypes.reindex(columns=group_l)
    if chains[0].startswith("TR"):
        return clonotypes, 'TCR'
    elif chains[0].startswith("IG"):
        return clonotypes, 'BCR'


class Split_tag(Step):
    """
    ## Features
    - Split scRNA-Seq fastq according to tag assignment.

    ## Output
    - `matrix/` Matrix files of each tag.(Optional)
    - `fastq/` Fastq files of each tag.(Optional)
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        if not (args.split_matrix or args.split_fastq or args.split_vdj or args.split_fl_vdj):
            return

        # set
        df_umi_tag = pd.read_csv(args.umi_tag_file, sep='\t', index_col=0)
        df_umi_tag = df_umi_tag.rename_axis('barcode').reset_index()
        self.tag_barcode_dict = {tag: set(row["barcode"].tolist()) for tag, row in df_umi_tag.groupby("tag")}

        if args.split_matrix:
            self.matrix_outdir = f'{args.outdir}/matrix/'
            if utils.check_arg_not_none(args, 'match_dir'):
                matrix_dir = utils.get_matrix_dir_from_match_dir(args.match_dir)
            elif utils.check_arg_not_none(args, 'matrix_dir'):
                matrix_dir = args.matrix_dir
            else:
                raise ValueError("--match_dir or --matrix_dir is required.")
            self.count_matrix = CountMatrix.from_matrix_dir(matrix_dir)

        if args.split_fastq:
            self.rna_fq_file = glob.glob(f'{args.match_dir}/*barcode/*_2.fq*')[0]

            fastq_outdir = f'{args.outdir}/fastqs/'
            os.system(f'mkdir -p {fastq_outdir}')

            self.r2_fastq_files_handle = {}
            self.r1_fastq_files_handle = {}
            for tag in self.tag_barcode_dict:
                r2_fastq_file_name = f'{fastq_outdir}/{tag}_2.fq'
                self.r2_fastq_files_handle[tag] = open(r2_fastq_file_name, 'w')
                r1_fastq_file_name = f'{fastq_outdir}/{tag}_1.fq'
                self.r1_fastq_files_handle[tag] = open(r1_fastq_file_name, 'w')

            self.tag_read_index_dict = defaultdict(set)

        if args.split_vdj:
            self.cell_confident_vdj = glob.glob(f'{args.vdj_dir}/*count_vdj/*cell_confident.tsv*')[0]

            self.vdj_outdir = f'{args.outdir}/vdj/'
            if not os.path.exists(self.vdj_outdir):
                os.system(f'mkdir -p {self.vdj_outdir}')
        
        if args.split_fl_vdj:
            self.fl_vdj_outdir = f'{args.outdir}/fl_vdj/'
            if not os.path.exists(self.fl_vdj_outdir):
                os.system(f'mkdir -p {self.fl_vdj_outdir}')
            
            # flv_CR
            if os.path.exists(f'{args.vdj_dir}/02.convert'):
                if os.path.exists(f'{args.vdj_dir}/03.assemble/match'):
                    match_out_dir = f'{args.vdj_dir}/03.assemble/match'
                else:
                    match_out_dir = f'{args.vdj_dir}/05.match'
                
                try: # old version
                    self.anno_file = glob.glob(f'{match_out_dir}/match_contigs.csv')[0]
                    self.fasta_file = glob.glob(f'{match_out_dir}/match_contig.fasta')[0]
                except IndexError: # latest version
                    self.anno_file = glob.glob(f'{match_out_dir}/matched_contig_annotations.csv')[0]
                    self.fasta_file = glob.glob(f'{match_out_dir}/matched_contig.fasta')[0]

            # flv_trust4
            else:
                self.anno_file = glob.glob(f'{args.vdj_dir}/04.summarize/*_filtered_contig.csv')[0]
                self.fasta_file = glob.glob(f'{args.vdj_dir}/04.summarize/*_filtered_contig.fasta')[0]


    @utils.add_log
    def write_r2_fastq_files(self):
        read_num = 0
        with pysam.FastxFile(self.rna_fq_file, 'r') as rna_fq:
            for read in rna_fq:
                read_num += 1
                attr = read.name.strip("@").split("_")
                barcode = attr[0]
                read_index = int(attr[2])
                for tag in self.tag_barcode_dict:
                    if barcode in self.tag_barcode_dict[tag]:
                        self.tag_read_index_dict[tag].add(read_index)
                        self.r2_fastq_files_handle[tag].write(str(read) + '\n')

                if read_num % 1000000 == 0:
                    self.write_r2_fastq_files.logger.info(f'{read_num} done')

        for tag in self.r2_fastq_files_handle:
            self.r2_fastq_files_handle[tag].close()

    @utils.add_log
    def write_r1_fastq_files(self):
        file_handles = [pysam.FastxFile(r1, 'r') for r1 in self.args.R1_read.split(',')]
        r1_read = itertools.chain(*file_handles)
        for read_index, read in enumerate(r1_read, start=1):
            for tag in self.tag_read_index_dict:
                if read_index in self.tag_read_index_dict[tag]:
                    self.r1_fastq_files_handle[tag].write(str(read) + '\n')

        for r1 in file_handles:
            r1.close()
        for tag in self.r1_fastq_files_handle:
            self.r1_fastq_files_handle[tag].close()

    @utils.add_log
    def split_matrix(self):
        for tag in self.tag_barcode_dict:
            tag_barcodes = list(self.tag_barcode_dict[tag])
            raw_barcodes = self.count_matrix.get_barcodes()
            tag_barcodes_indices = [raw_barcodes.index(barcode) for barcode in tag_barcodes]
            tag_barcodes_indices.sort()
            slice_matrix = self.count_matrix.slice_matrix(tag_barcodes_indices)

            tag_matrix_dir = f'{self.matrix_outdir}/{tag}_{FILTERED_MATRIX_DIR_SUFFIX[0]}/'
            slice_matrix.to_matrix_dir(tag_matrix_dir)

    @utils.add_log
    def split_vdj(self):

        df_vdj = pd.read_csv(self.cell_confident_vdj, sep='\t')
        for tag in self.tag_barcode_dict:
            tag_barcodes = list(self.tag_barcode_dict[tag])
            temp = df_vdj[df_vdj.barcode.isin(tag_barcodes)]
            if not temp.empty:
                clonotypes, seqtype = get_clonotypes_table(temp)
                temp.to_csv(f'{self.vdj_outdir}/{tag}_{seqtype}_cell_confident.tsv', sep='\t', index=False)
                clonotypes.to_csv(f'{self.vdj_outdir}/{tag}_{seqtype}_clonotypes.tsv', sep='\t', index=False)
            else:
                continue
    
    @utils.add_log
    def split_fl_vdj(self):
        df_anno = pd.read_csv(self.anno_file, keep_default_na=False)
        if df_anno.chain[0].startswith('IG'):
            seqtype = 'BCR'
        else:
            seqtype = 'TCR'

        for tag in self.tag_barcode_dict:
            tag_barcodes = set(self.tag_barcode_dict[tag])
            df_temp = df_anno[df_anno.barcode.isin(tag_barcodes)]
            if not df_temp.empty:
                df_temp.to_csv(f'{self.fl_vdj_outdir}/{tag}_{seqtype}_contig_annotations.csv', sep=',', index=False)
                fasta_temp = open(f'{self.fl_vdj_outdir}/{tag}_{seqtype}_contig.fasta' ,'w')
                with pysam.FastxFile(self.fasta_file) as raw_fasta:
                    for entry in raw_fasta:
                        name = entry.name
                        attrs = name.split('_')
                        cb = attrs[0]
                        if cb in tag_barcodes:
                            new_name = cb + '_' + attrs[1] + '_' + attrs[2]
                            seq = entry.sequence
                            fasta_temp.write(f'>{new_name}\n{seq}\n')
                fasta_temp.close()

                if 'clonotype_id' in df_temp.columns:
                    df_temp = df_temp[df_temp.clonotype_id != ''] # del clonotype for trust4
                split_clonotypes = f'{self.fl_vdj_outdir}/{tag}_{seqtype}_clonotypes.csv'
                self.get_fl_clonotypes_table(df_temp, split_clonotypes)

    @staticmethod
    def get_fl_clonotypes_table(df_temp, split_clonotypes):
        df_match = df_temp[df_temp['productive'] == True]
        df_match['chain_cdr3aa'] = df_match[['chain', 'cdr3']].apply(':'.join, axis=1)

        match_clonotypes = open(split_clonotypes, 'w')
        match_clonotypes.write('barcode\tcdr3s_aa\n')
        for cb in set(df_match.barcode):
            temp = df_match[df_match['barcode']==cb].sort_values(by='chain', ascending=True)
            chain_pair = ';'.join(temp['chain_cdr3aa'].tolist())
            match_clonotypes.write(f'{cb}\t{chain_pair}\n')
        match_clonotypes.close()

        df_match_clonetypes = pd.read_csv(split_clonotypes, sep='\t', index_col=None)
        df_match_clonetypes = df_match_clonetypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
        df_match_clonetypes.rename(columns={'barcode': 'frequency'}, inplace=True)
        sum_f = df_match_clonetypes['frequency'].sum()
        df_match_clonetypes['proportion'] = df_match_clonetypes['frequency'].apply(lambda x: x/sum_f)
        df_match_clonetypes.sort_values(by='frequency', ascending=False, inplace=True)
        df_match_clonetypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_match_clonetypes.shape[0]+1)]
        df_match_clonetypes = df_match_clonetypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
        df_match_clonetypes.to_csv(split_clonotypes, sep=',', index=False)

    @utils.add_log
    def run(self):
        if self.args.split_matrix:
            self.split_matrix()
        if self.args.split_fastq:
            self.write_r2_fastq_files()
            self.write_r1_fastq_files()
        if self.args.split_vdj:
            self.split_vdj()
        if self.args.split_fl_vdj:
            self.split_fl_vdj()


def split_tag(args):
    with Split_tag(args) as runner:
        runner.run()


def get_opts_split_tag(parser, sub_program):
    parser.add_argument(
        "--split_fastq",
        help="If used, will split scRNA-Seq fastq file according to tag assignment.",
        action='store_true',
    )
    parser.add_argument(
        "--split_matrix",
        help="If used, will split scRNA-Seq matrix file according to tag assignment.",
        action='store_true',
    )
    parser.add_argument(
        "--split_vdj",
        help="If used, will split scRNA-Seq vdj count file according to tag assignment.",
        action='store_true',
    )
    parser.add_argument(
        "--split_fl_vdj",
        help='If used, will split scRNA-Seq full-length vdj annotation, fasta, clonotypes file according to tag assignment.',
        action='store_true',
    )
    parser.add_argument("--vdj_dir", help="Match celescope vdj directory. Required when --split_vdj or --split_fl_vdj is specified.")
    if sub_program:
        parser.add_argument("--umi_tag_file", help="UMI tag file.", required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--matrix_dir", help="Match celescope scRNA-Seq matrix directory.")
        parser.add_argument("--R1_read", help='R1 read path.')
        s_common(parser)
