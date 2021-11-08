import subprocess
from collections import defaultdict
from multiprocessing import Pool
from itertools import groupby
from functools import partial
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
import numpy as np
import pyranges as pr
import pysam
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

import celescope.tools.utils as utils
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.rna.mkref import parse_genomeDir_rna


OTSU_READ_MIN = 30



def parse_vcf(vcf_file, cols=('chrom', 'pos', 'alleles',), infos=('VID',)):
    '''
    parse vcf into df
    '''
    vcf = pysam.VariantFile(vcf_file)
    rec_dict_list = []
    for rec in vcf.fetch():
        rec_dict = {}
        for col in cols:
            rec_dict[col] = getattr(rec, col)
            # if ref == alt: alleles=(ref,)
            # else alleles=(ref, alt)
            if col == 'alleles':
                rec_dict['ref'] = rec_dict['alleles'][0]
                rec_dict['alt'] = '.'
                if len(rec_dict['alleles']) >= 2:
                    rec_dict['alt'] = ','.join(rec_dict['alleles'][1:])

        for info in infos:
            rec_dict[info] = rec.info[info]
        rec_dict_list.append(rec_dict)

    df = pd.DataFrame.from_dict(rec_dict_list)
    vcf.close()
    return df


def read_CID(CID_file):
    df_index = pd.read_csv(CID_file, sep='\t', index_col=0).reset_index()
    df_valid = df_index[df_index['valid'] == True]
    return df_index, df_valid




def map_vcf_row(row, df_cell_vcf):
    """
    get ref and UMI for each variant
    Args:
        row: each row from merged_vcf
    """
    pos = row['pos']
    chrom = row['chrom']
    alt = row['alt']
    df_pos = df_cell_vcf[(df_cell_vcf['pos'] == pos) & (df_cell_vcf['chrom'] == chrom)]
    df_ref = df_pos[df_pos['alt'] == '.']
    df_alt = df_pos[df_pos['alt'] == alt]
    ref_UMI = 0
    alt_UMI = 0
    if df_ref.shape[0] != 0:
        ref_UMI = get_DP4(df_ref, 'ref')
    if df_alt.shape[0] != 0:
        alt_UMI = get_DP4(df_alt, 'alt')
    return ref_UMI, alt_UMI


def get_DP4(row, alt):
    DP4 = row['DP4'].iloc[0]
    if alt == 'ref':
        indexs = [0, 1]
    elif alt == 'alt':
        indexs = [2, 3]
    umi = sum([DP4[index] for index in indexs])
    return umi


class Variant_calling(Step):
    """
    Features
    - Perform variant calling at single cell level.

    Output

    `{sample}_VID.tsv` A unique numeric ID is assigned for each variant, 

    - `RID`: Target region ID. This column will be added when `--panel` option were provided.

    `{sample}_CID.tsv` A unique numeric ID is assigned for each cell.

    `{sample}_RID.tsv` A unique numeric ID is assigned for each target region. This file will be created when `--panel` option were provided.


    `{sample}_merged.vcf ` VCF file containing all variants of all cells. `VID` and `CID` are added to the `INFO` column.

    `{sample}_filter.vcf` VCF file after filtering. Invalid `CID`s are removed from the `INFO` column.

    `{sample}_variant_count.tsv`  Reference and variant supporting reads/UMIs counts.

    `{sample}_filter_variant_count.tsv`  Reference and variant supporting reads/UMIs counts after filtering.

    `{sample}_support.mtx` Support matrix in [Matrix Market Exchange Formats](https://math.nist.gov/MatrixMarket/formats.html). Rows 
    are variants(VID) and columns are cells(CID). The value can be 1, 2 or 3.
    
    1 : all reads/UMIs at the position support the ref allele.  
    2 : all reads/UMIs at the position support the alt allele.  
    3 : one or more reads/UMIs support both the alt and the ref allele.  
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # set
        self.barcodes, _num = utils.read_barcode_file(args.match_dir)
        self.fasta = parse_genomeDir_rna(args.genomeDir)['fasta']
        self.df_vcf = None
        self.panel = args.panel
        self.bed = utils.get_bed_file_path(self.panel)

        # out
        self.splitN_bam = f'{self.out_prefix}_splitN.bam'
        self.splitN_bam_name_sorted = f'{self.out_prefix}_splitN_name_sorted.bam'
        self.CID_file = f'{self.out_prefix}_CID.tsv'
        self.VID_file = f'{self.out_prefix}_VID.tsv'
        self.RID_file = f'{self.out_prefix}_RID.tsv'
        self.cells_dir = f'{self.outdir}/cells/'
        self.merged_vcf_file = f'{self.out_prefix}_merged.vcf'
        self.filter_vcf_file = f'{self.out_prefix}_filter.vcf'
        self.variant_count_file = f'{self.out_prefix}_variant_count.tsv'
        self.otsu_dir = f'{self.out_prefix}_otsu/'
        self.otsu_threshold_file = f'{self.otsu_dir}/{self.sample}_otsu_threshold.tsv'
        self.filter_variant_count_file = f'{self.out_prefix}_filter_variant_count.tsv'
        if args.min_support_read == 'auto':
            utils.check_mkdir(self.otsu_dir)
        self.support_matrix_file = f'{self.out_prefix}_support.mtx'

        self.raw_bcf_file = f'{self.out_prefix}_raw.bcf'
        self.raw_vcf_file = f'{self.out_prefix}_raw.vcf'

        # run
        self.CID_dict = defaultdict(dict)
        for index, barcode in enumerate(self.barcodes):
            CID = index + 1
            self.CID_dict[barcode]['CID'] = CID
            self.CID_dict[barcode]['valid'] = False


    @utils.add_log
    def SplitNCigarReads(self):
        cmd = (
            f'gatk '
            f'SplitNCigarReads '
            f'--do-not-fix-overhangs '
            f'-R {self.fasta} '
            f'-I {self.args.bam} '
            f'-O {self.splitN_bam} '
        )
        Variant_calling.SplitNCigarReads.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def call_variants(self):
        cmd = (
            f'bcftools mpileup '
            f'-f {self.fasta} '
            f'--threads {self.thread} '
            f'--annotate DP,AD -d 10000 '
            f'-o {self.raw_bcf_file} '
            f'{self.splitN_bam} '
        )
        if self.bed:
            cmd += f' --regions-file {self.bed} '
        self.debug_subprocess_call(cmd)

        cmd = (
            f'bcftools call '
            f'-mv -Ov '
            f'-o {self.raw_vcf_file} '
            f'{self.raw_bcf_file} '
        )
        self.debug_subprocess_call(cmd)

    def write_CID_file(self):
        # out CID file
        df_CID = pd.DataFrame(self.CID_dict).T
        df_CID.index.name = 'barcode'
        df_CID.to_csv(self.CID_file, sep='\t')


    def read_CID(self):
        return read_CID(self.CID_file)

    @utils.add_log
    def merge_vcf(self):
        '''
        merge cell vcf into one non-duplicated vcf
        add VID(variant ID) and CID(cell ID)
        '''
        _df_index, df_valid = self.read_CID()
        CIDs = list(df_valid['CID'])
        # variant dict
        v_cols = ['chrom', 'pos', 'alleles']
        v_dict = dict()

        for CID in CIDs:
            CID = str(CID)
            vcf_file = f'{self.outdir}/cells/cell{CID}/cell{CID}_norm.vcf'
            with pysam.VariantFile(vcf_file, 'r') as vcf:
                for rec in vcf.fetch():
                    v = ','.join([str(getattr(rec, col)) for col in v_cols])
                    if not v in v_dict:
                        v_dict[v] = dict()
                        for col in v_cols:
                            v_dict[v][col] = getattr(rec, col)
                        v_dict[v]['CID'] = [CID]
                        v_dict[v]['record'] = rec
                    else:
                        v_dict[v]['CID'].append(CID)

        # output
        def get_vcf_header(CIDs):
            CID = CIDs[0]
            vcf_file = f'{self.outdir}/cells/cell{CID}/cell{CID}_norm.vcf'
            with pysam.VariantFile(vcf_file, 'r') as vcf:
                return vcf.header

        vcf_header = get_vcf_header(CIDs)
        vcf_header.info.add('VID', number=1, type='String', description='Variant ID')
        vcf_header.info.add('CID', number=1, type='String', description='Cell ID')
        with pysam.VariantFile(self.merged_vcf_file, 'w', header=vcf_header) as merged_vcf:
            VID = 0
            for _keys, values in sorted(v_dict.items(), 
                key=lambda x: (x[1]['chrom'], x[1]['pos'], x[1]['alleles'])):
                VID += 1
                rec = values['record']
                CID = ','.join(values['CID'])
                record = merged_vcf.new_record()
                cols = ['chrom', 'pos', 'alleles']
                for col in cols:
                    setattr(record, col, getattr(rec, col))
                record.info['VID'] = str(VID)
                record.info['CID'] = CID
                merged_vcf.write(record)

    @utils.add_log
    def write_VID_file(self):
        df_vcf = parse_vcf(self.merged_vcf_file)
        df_VID = df_vcf.loc[:, ['VID', 'chrom', 'pos', 'ref', 'alt']]

        df_VID.to_csv(self.VID_file, sep='\t', index=False)

    def otsu_threshold(self, array):
        threshold = utils.otsu_min_support_read(array, self.otsu_plot)
        return threshold
   
    def get_bed_df(self):
        bed_df = utils.get_gene_region_from_bed(self.panel)[1]
        return bed_df

    @utils.add_log
    def write_RID_file(self):
        rid_file = self.get_bed_df()
        rid_file.insert(0, 'RID', [rid + 1 for rid in range(len(rid_file))], allow_duplicates=False)
        rid_file.to_csv(self.RID_file,sep = '\t',index = False)
    
    @utils.add_log
    def get_region_dict(self):
        """
        Return a dict
        - key: vid
        - value: rid
        """
        df_rid = pd.read_table(self.RID_file,sep = '\t')
        df_vid = pd.read_table(self.VID_file,sep = '\t')
        
        bed_file_df = self.get_bed_df()
        gr = pr.PyRanges(bed_file_df)
        region_dict = {}
        for _, row in df_vid.iterrows():
            vid = row['VID']
            bed_region = gr[str(row['chrom']), row['pos']: row['pos'] + 1]
            if bed_region:
                bed_chr = bed_region.as_df().loc[:,"Chromosome"].astype(int).to_list()
                bed_start = bed_region.as_df().loc[:,"Start"].to_list()
                bed_end = bed_region.as_df().loc[:,"End"].to_list()
                #get rid
                for chrs, start, end in zip(bed_chr,bed_start,bed_end):
                    rid = df_rid[(df_rid.loc[:,"Chromosome"] == chrs) 
                                 & (df_rid.loc[:,"Start"] == start) 
                                 & (df_rid.loc[:,"End"] == end)].loc[:,"RID"].to_list()[0]
                    region_dict[vid] = rid
            else:
                region_dict[vid] = "None"
        return region_dict
    
    @utils.add_log
    def add_RID(self):
        df_vid = pd.read_table(self.VID_file,sep = '\t')
        df_filter_variant_conut = pd.read_table(self.filter_variant_count_file,sep = '\t')
        
        region_dict = self.get_region_dict()
        df_vid['RID'] = df_vid['VID'].apply(lambda x: region_dict[x])
        df_filter_variant_conut['RID'] = df_filter_variant_conut['VID'].apply(lambda x: region_dict[x])

        df_vid.to_csv(self.VID_file,sep = '\t',index = None)
        df_filter_variant_conut.to_csv(self.filter_variant_count_file,sep = '\t',index = None)
    
    @utils.add_log
    def filter_vcf(self):
        """
        filter cells with zero variant UMI
        """
        df_filter = pd.read_csv(self.filter_variant_count_file, sep='\t')
        df_filter = df_filter[df_filter['alt_count']>0]

        vcf = pysam.VariantFile(self.merged_vcf_file, 'r')
        vcf_header = vcf.header
        filter_vcf = pysam.VariantFile(self.filter_vcf_file, 'w', header=vcf_header)
        for rec in vcf.fetch():
            VID = int(rec.info['VID'])
            CIDs = df_filter[df_filter['VID']==VID]['CID'].values
            rec.info['CID'] = ','.join([str(CID) for CID in CIDs])
            if rec.info['CID']:
                filter_vcf.write(rec)
        filter_vcf.close()
        vcf.close()


    @utils.add_log
    def write_support_matrix(self):
        def set_support_bit(row):
            ref_bit = 1 if row['ref_count'] > 0 else 0
            alt_bit = 2 if row['alt_count'] > 0 else 0
            support_bit = ref_bit + alt_bit
            return support_bit

        df_variant_count = pd.read_csv(self.filter_variant_count_file, sep='\t')
        df_variant_count['support'] = df_variant_count.apply(set_support_bit, axis=1)
        support_mtx = coo_matrix(
            (df_variant_count.support, (df_variant_count.VID - 1, df_variant_count.CID - 1))
        )
        mmwrite(self.support_matrix_file, support_mtx)

    def run(self):

        self.SplitNCigarReads()
        self.call_variants()
        #self.clean_up()


@utils.add_log
def variant_calling(args):

    step_name = 'variant_calling'
    variant_calling_obj = Variant_calling(args, step_name)
    variant_calling_obj.run()


def get_opts_variant_calling(parser, sub_program):

    parser.add_argument("--genomeDir", help=HELP_DICT['genomeDir'], required=True)
    parser.add_argument(
        "--min_support_read",
        help="""Minimum number of reads support a variant. If `auto`(default), otsu method will be used to determine this value.""",
        default='auto',
    )
    parser.add_argument(
        "--panel",
        help = "The prefix of bed file, such as `lung_1`.",
        default = ''
        )
    if sub_program:
        parser.add_argument(
            "--bam",
            help='Input BAM file from step `target_metrics`. ',
            required=True
        )
        parser.add_argument(
            "--match_dir",
            help=HELP_DICT['match_dir'],
            required=True
        )
        s_common(parser)
