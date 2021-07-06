import logging
import os
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
import pysam
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

import celescope.tools.utils as utils
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.rna.mkref import parse_genomeDir_rna


def parse_vcf(vcf_file, cols=('chrom', 'pos', 'alleles',), infos=('VID',)):
    '''
    parse vcf into df
    '''
    vcf = pysam.VariantFile(vcf_file)
    df = pd.DataFrame(columns=list(cols) + list(infos))
    rec_dict = {}
    for rec in vcf.fetch():

        for col in cols:
            rec_dict[col] = getattr(rec, col)
            # if ref == alt: alleles=(ref,)
            # else alleles=(ref, alt)
            if col == 'alleles':
                rec_dict['ref'] = rec_dict['alleles'][0]
                rec_dict['alt'] = '.'
                if len(rec_dict['alleles']) == 2:
                    rec_dict['alt'] = rec_dict['alleles'][1]

        for info in infos:
            rec_dict[info] = rec.info[info]

        df = df.append(pd.Series(rec_dict), ignore_index=True)
    return df


def read_CID(CID_file):
    df_index = pd.read_csv(CID_file, sep='\t', index_col=0, dtype=object)
    df_valid = df_index[df_index['valid'] == 'True']
    return df_index, df_valid


class Variant_calling(Step):
    """
    Features
    - Perform variant calling.

    Output

    `{sample}_VID.tsv` A unique numeric ID is assigned for each variant.

    `{sample}_CID.tsv` A unique numeric ID is assigned for each cell.

    `{sample}_variant_count.tsv`  Reference and variant supporting reads/UMIs count.

    `{sample}_support.mtx` Support matrix, only high quality bases are considered.   
    0 : no reads/UMIs cover the position.  
    1 : all reads/UMIs at the position support the ref allele.  
    2 : all reads/UMIs at the position support the alt allele.  
    3 : one or more reads/UMIs support both the alt and the ref allele.  
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # set
        self.barcodes, _num = utils.read_barcode_file(args.match_dir)
        self.fasta = parse_genomeDir_rna(args.genomeDir)['fasta']
        if args.vcf:
            self.vcf_bool = True
        else:
            self.vcf_bool = False
        self.df_vcf = None

        # out
        self.splitN_bam = f'{self.out_prefix}_splitN.bam'
        self.CID_file = f'{self.out_prefix}_CID.tsv'
        self.VID_file = f'{self.out_prefix}_VID.tsv'
        self.final_vcf_file = f'{self.out_prefix}.vcf'
        self.variant_count_file = f'{self.out_prefix}_variant_count.tsv'
        self.support_matrix_file = f'{self.out_prefix}_support.mtx'

    @utils.add_log
    def SplitNCigarReads(self):
        cmd = (
            f'gatk '
            f'SplitNCigarReads '
            f'-R {self.fasta} '
            f'-I {self.args.bam} '
            f'-O {self.splitN_bam} '
        )
        Variant_calling.SplitNCigarReads.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def split_bam(self):
        '''
        input:
            bam: bam from splitN
            barcodes: cell barcodes, list
        ouput:
            bam_dict: assign reads to cell barcodes and UMI
            count_dict: UMI counts per cell
            CID: assign ID(1-based) to cells
        '''

        # init
        bam_dict = defaultdict(list)
        CID_dict = defaultdict(dict)
        cells_dir = f'{self.outdir}/cells/'

        # read bam and split
        samfile = pysam.AlignmentFile(self.splitN_bam, "rb")
        header = samfile.header
        for read in samfile:
            try:
                barcode = read.get_tag('CB')
            except KeyError:
                continue
            if barcode in self.barcodes:
                CID = self.barcodes.index(barcode) + 1
                read.set_tag(tag='CL', value=f'CELL{CID}', value_type='Z')

                # assign read to barcode
                bam_dict[barcode].append(read)

        self.split_bam.logger.info('writing cell bam...')
        # write new bam
        CID = 0
        for barcode in self.barcodes:
            # init
            CID += 1
            CID_dict[CID]['barcode'] = barcode
            CID_dict[CID]['valid'] = False

            # out bam
            if barcode in bam_dict:
                cell_dir = f'{cells_dir}/cell{CID}'
                cell_bam_file = f'{cell_dir}/cell{CID}.bam'
                if not os.path.exists(cell_dir):
                    os.makedirs(cell_dir)
                CID_dict[CID]['valid'] = True
                cell_bam = pysam.AlignmentFile(
                    f'{cell_bam_file}', "wb", header=header)
                for read in bam_dict[barcode]:
                    cell_bam.write(read)
                cell_bam.close()

        # out CID
        df_CID = pd.DataFrame(CID_dict).T
        df_CID.index.name = 'CID'
        df_CID.to_csv(self.CID_file, sep='\t')

    @staticmethod
    @utils.add_log
    def call_snp(CID, outdir, fasta):

        Variant_calling.call_snp.logger.info('Processing Cell %s' % CID)
        bam = f'{outdir}/cells/cell{CID}/cell{CID}.bam'
        # sort
        sorted_bam = f'{outdir}/cells/cell{CID}/cell{CID}_sorted.bam'
        cmd_sort = (
            f'samtools sort {bam} -o {sorted_bam}'
        )
        subprocess.check_call(cmd_sort, shell=True)

        # mpileup
        bcf = f'{outdir}/cells/cell{CID}/cell{CID}.bcf'
        cmd_mpileup = (
            f'bcftools mpileup -Ou '
            f'-f {fasta} '
            f'{sorted_bam} -o {bcf} '
        )
        subprocess.check_call(cmd_mpileup, shell=True)

        # call
        out_vcf = f'{outdir}/cells/cell{CID}/cell{CID}.vcf'
        cmd_call = (
            f'bcftools call -mv -Ov '
            f'-o {out_vcf} '
            f'{bcf}'
            f'>/dev/null 2>&1 '
        )
        subprocess.check_call(cmd_call, shell=True)

        # norm
        norm_vcf = f'{outdir}/cells/cell{CID}/cell{CID}_norm.vcf'
        cmd_norm = (
            f'bcftools norm -d none '
            f'-f {fasta} '
            f'{out_vcf} '
            f'-o {norm_vcf} '
        )
        subprocess.check_call(cmd_norm, shell=True)

        # call all position
        out_all_vcf = f'{outdir}/cells/cell{CID}/cell{CID}_all.vcf'
        cmd_all_call = (
            f'bcftools call -m -Ov '
            f'-o {out_all_vcf} '
            f'{bcf}'
            f'>/dev/null 2>&1 '
        )
        subprocess.check_call(cmd_all_call, shell=True)

        # norm all
        norm_all_vcf = f'{outdir}/cells/cell{CID}/cell{CID}_all_norm.vcf'
        cmd_all_norm = (
            f'bcftools norm -d none '
            f'-f {fasta} '
            f'{out_all_vcf} '
            f'-o {norm_all_vcf} '
        )
        subprocess.check_call(cmd_all_norm, shell=True)

    @utils.add_log
    def call_all_snp(self):
        all_res = []
        _df_index, df_valid = self.read_CID()
        CID_arg = df_valid.index
        outdir_arg = [self.outdir] * len(CID_arg)
        fasta_arg = [self.fasta] * len(CID_arg)
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(self.call_snp, CID_arg, outdir_arg, fasta_arg):
                all_res.append(res)

    def read_CID(self):
        return read_CID(self.CID_file)

    @utils.add_log
    def merge_vcf(self):
        '''
        if vcf not provided,
        merge cell vcf into one non-duplicated vcf
        add VID(variant ID) and CID(cell ID)
        '''
        _df_index, df_valid = self.read_CID()
        CIDs = df_valid.index

        # variant dict
        v_cols = ['chrom', 'pos', 'alleles']
        v_dict = {}

        for CID in CIDs:
            CID = str(CID)
            vcf_file = f'{self.outdir}/cells/cell{CID}/cell{CID}_norm.vcf'
            vcf = pysam.VariantFile(vcf_file, 'r')
            for rec in vcf.fetch():
                v = ','.join([str(getattr(rec, col)) for col in v_cols])
                if not v in v_dict:
                    v_dict[v] = dict()
                    v_dict[v]['CID'] = [CID]
                    v_dict[v]['record'] = rec
                else:
                    v_dict[v]['CID'].append(CID)

        # output
        def get_vcf_header(CIDs):
            CID = CIDs[0]
            vcf_file = f'{self.outdir}/cells/cell{CID}/cell{CID}_norm.vcf'
            vcf = pysam.VariantFile(vcf_file, 'r')
            return vcf.header
        vcf_header = get_vcf_header(CIDs)
        vcf_header.info.add('VID', number=1, type='String', description='Variant ID')
        vcf_header.info.add('CID', number=1, type='String', description='Cell ID')
        merged_vcf = pysam.VariantFile(self.final_vcf_file, 'w', header=vcf_header)

        VID = 0
        for v in sorted(v_dict.keys()):
            VID += 1
            rec = v_dict[v]['record']
            CID = ','.join(v_dict[v]['CID'])
            record = merged_vcf.new_record()
            cols = ['chrom', 'pos', 'alleles']
            for col in cols:
                setattr(record, col, getattr(rec, col))
            record.info['VID'] = str(VID)
            record.info['CID'] = CID
            merged_vcf.write(record)
        merged_vcf.close()

    @utils.add_log
    def write_VID_file(self):
        df_vcf = parse_vcf(self.final_vcf_file)
        df_VID = df_vcf.loc[:, ['VID', 'chrom', 'pos', 'ref', 'alt']]
        df_VID.to_csv(self.VID_file, sep='\t', index=False)

    @utils.add_log
    def add_VID(self):
        vcf = pysam.VariantFile(self.args.vcf, 'r')
        vcf_header = vcf.header
        if 'VID' in vcf_header.info:
            logging.info('VID is already in vcf file!')
            return
        vcf_header.info.add('VID', number=1, type='String', description='Variant ID')
        VID_vcf = pysam.VariantFile(self.final_vcf_file, 'w', header=vcf_header)
        VID = 0
        for rec in vcf.fetch():
            VID += 1
            rec.info['VID'] = str(VID)
            VID_vcf.write(rec)
        VID_vcf.close()

    @staticmethod
    def cell_UMI(CID, outdir, final_vcf_file):
        df_vcf = parse_vcf(final_vcf_file)
        df_UMI = pd.DataFrame(columns=['VID', 'CID', 'ref_count', 'alt_count'])
        norm_all_vcf = f'{outdir}/cells/cell{CID}/cell{CID}_all_norm.vcf'
        df_cell_vcf = parse_vcf(norm_all_vcf, infos=['DP4'])

        def get_DP4(row, alt):
            DP4 = row['DP4'].iloc[0]
            if alt == 'ref':
                indexs = [0, 1]
            elif alt == 'alt':
                indexs = [2, 3]
            umi = sum([DP4[index] for index in indexs])
            return umi

        def map_vcf_row(row, df_cell_vcf):
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
            return ref_UMI, alt_UMI, pos, chrom, alt

        for index in df_vcf.index:
            row = df_vcf.loc[index, ]
            ref_UMI, alt_UMI, _pos, _chrom, _alt = map_vcf_row(row, df_cell_vcf)
            if (ref_UMI + alt_UMI) != 0:
                VID = row['VID']
                dic = {
                    'VID': VID,
                    'CID': CID,
                    'ref_count': ref_UMI,
                    'alt_count': alt_UMI,
                }
                df_UMI = df_UMI.append(dic, ignore_index=True)
        return df_UMI

    @utils.add_log
    def get_UMI(self):
        '''
        get variant and ref UMI supporting an allele
        '''
        _df_index, df_valid = self.read_CID()

        df_UMI_list = []
        CID_arg = list(df_valid.index)
        outdir_arg = [self.outdir] * len(CID_arg)
        final_vcf_file_arg = [self.final_vcf_file] * len(CID_arg)
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(Variant_calling.cell_UMI, CID_arg, outdir_arg, final_vcf_file_arg):
                df_UMI_list.append(res)

        df_UMI = pd.concat(df_UMI_list)
        df_UMI['VID'] = df_UMI['VID'].astype('int')
        df_UMI.sort_values(by=['VID', 'CID'], inplace=True)
        df_UMI.to_csv(self.variant_count_file, sep='\t', index=False)

    @utils.add_log
    def write_support_matrix(self):
        def set_support_bit(row):
            ref_bit = 1 if row['ref_count'] > 0 else 0
            alt_bit = 2 if row['alt_count'] > 0 else 0
            support_bit = ref_bit + alt_bit
            return support_bit

        df_variant_count = pd.read_csv(self.variant_count_file, sep='\t')
        df_variant_count['support'] = df_variant_count.apply(set_support_bit, axis=1)
        support_mtx = coo_matrix(
            (df_variant_count.support, (df_variant_count.VID - 1, df_variant_count.CID - 1))
        )
        mmwrite(self.support_matrix_file, support_mtx)

    def run(self):
        self.SplitNCigarReads()
        self.split_bam()
        self.call_all_snp()
        if self.vcf_bool:
            self.add_VID()
        else:
            self.merge_vcf()
        self.write_VID_file()
        self.get_UMI()
        self.write_support_matrix()
        self.clean_up()


@utils.add_log
def variant_calling(args):

    step_name = 'variant_calling'
    variant_calling_obj = Variant_calling(args, step_name)
    variant_calling_obj.run()


def get_opts_variant_calling(parser, sub_program):

    parser.add_argument("--genomeDir", help=HELP_DICT['genomeDir'], required=True)
    parser.add_argument(
        "--vcf",
        help="""VCF file. If vcf file is not provided, celescope will perform variant calling at single cell level 
and use these variants as input vcf.""",
        required=False
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
