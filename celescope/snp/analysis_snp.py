
import configparser
import subprocess

import pandas as pd
import pysam

import celescope.tools.utils as utils
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common
from celescope.snp.variant_calling import read_CID
from celescope.__init__ import HELP_DICT


class Analysis_variant(Step, AnalysisMixin):

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        AnalysisMixin.__init__(self, args)
        self.filter_variant_count_file = args.filter_variant_count_file
        self.CID_file = args.CID_file
        self.filter_vcf = args.filter_vcf
        self.annovar_config = args.annovar_config
        self.match_dir = args.match_dir
        self.vcf_GT = None

    def get_df_count_tsne(self):
        '''
        output: f'{self.outdir}/{self.sample}_count_tsne.tsv'
        '''
        df_vc = pd.read_csv(self.filter_variant_count_file, sep='\t')
        df_vc = df_vc[df_vc["alt_count"] > 0]
        df_vc_cell = df_vc.groupby('CID').agg({
            'alt_count': 'count',
            'VID': list,
        })

        df_CID, _df_valid = read_CID(self.CID_file)
        df_CID = df_CID.reset_index()
        tsne_df_CID = pd.merge(self.tsne_df, df_CID, on='barcode', how='left')

        df_vc_barcode = pd.merge(df_vc_cell, df_CID, on='CID')
        df_vc_barcode_tsne = pd.merge(df_vc_barcode, tsne_df_CID, on=['barcode', 'CID'], how='right')
        df_vc_barcode_tsne['value'] = df_vc_barcode_tsne['alt_count']
        df_vc_barcode_tsne['value'] = df_vc_barcode_tsne['value'].fillna(0)
        df_vc_barcode_tsne['value'].astype('int32')
        df_count_tsne = df_vc_barcode_tsne
        # out
        out_file = f'{self.outdir}/{self.sample}_count_tsne.tsv'
        df_out = df_vc_barcode_tsne[df_vc_barcode_tsne['value'] > 0]
        cols = ['CID', 'alt_count', 'VID', 'barcode', 'cluster', 'tSNE_1', 'tSNE_2', 'Gene_Counts']
        df_out = df_out[cols]
        df_out.to_csv(out_file, sep='\t')
        return df_count_tsne

    def get_count_tsne(self, df_count_tsne):
        def return_text(row):
            text = f'CID:{row["CID"]} <br>Variants:{str(int(row["value"]))} <br>VID:{row["VID"]}'
            return text
        tSNE_1 = list(df_count_tsne.tSNE_1)
        tSNE_2 = list(df_count_tsne.tSNE_2)
        text = list(df_count_tsne.apply(return_text, axis=1))
        value = list(df_count_tsne.value)
        title = 't-SNE plot Colored by Cell Variant Counts'
        count_tsne = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "text": text, 'value': value, 'title': title}
        return count_tsne

    def add_GT(self):
        '''
        add genotype to VCF file to avoid vcf parse error
        '''
        vcf = pysam.VariantFile(self.filter_vcf, 'r')
        out_vcf_file = f'{self.outdir}/{self.sample}_addGT.vcf'
        out_vcf = pysam.VariantFile(out_vcf_file, 'w', header=vcf.header)
        for rec in vcf:
            for sample in rec.samples:
                rec.samples[sample]["GT"] = (1, 1)
                out_vcf.write(rec)
        vcf.close()
        out_vcf.close()
        self.vcf_GT = out_vcf_file

    def get_df_table(self):

        df_vcf = utils.parse_vcf(self.vcf_GT, infos=['VID', 'CID'])
        df_annovar = self.annovar()
        df_vcf = pd.concat((df_vcf, df_annovar), axis=1)
        df_vcf["nCell"] = df_vcf["CID"].apply(func=lambda row: 1 if isinstance(row, str) else len(row))

        out_df_vcf = f'{self.outdir}/{self.sample}_variant_table.tsv'
        df_vcf.to_csv(out_df_vcf, sep='\t', index=False)

        cols = ['VID', 'Chrom', 'Pos', 'Alleles', 'Gene', 'nCell', 'mRNA', 'Protein', 'COSMIC']
        df_vcf = df_vcf[cols]
        return df_vcf

    def run(self):
        self.add_GT()
        cluster_tsne = self.get_cluster_tsne(colname='cluster', tsne_df=self.tsne_df)
        df_count_tsne = self.get_df_count_tsne()
        count_tsne = self.get_count_tsne(df_count_tsne)
        df_vcf = self.get_df_table()
        table_dict = Step.get_table(title='Variant table', table_id='variant_table', df_table=df_vcf)

        self.add_data_item(cluster_tsne=cluster_tsne)
        self.add_data_item(count_tsne=count_tsne)
        self.add_data_item(table_dict=table_dict)
        self.clean_up()

    @utils.add_log
    def annovar(self):

        # config
        config = configparser.ConfigParser()
        config.read(self.annovar_config)
        section = config['ANNOVAR']
        annovar_dir = section['dir']
        db = section['db']
        buildver = section['buildver']
        protocol = section['protocol']
        operation = section['operation']

        # convert
        input_file = f'{self.outdir}/{self.sample}.input'
        cmd = (
            f'perl {annovar_dir}/convert2annovar.pl '
            f'-format vcf4 '
            f'--includeinfo '
            f'{self.vcf_GT} > {input_file}'
        )
        subprocess.check_call(cmd, shell=True)

        # annotate
        cmd = (
            f'perl {annovar_dir}/table_annovar.pl '
            f'{input_file} '
            f'{db} '
            f'-buildver {buildver} '
            f'-protocol {protocol} '
            f'-operation {operation} '
            f'-out {self.outdir}/{self.sample} '
            f'--otherinfo '
        )
        Analysis_variant.annovar.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

        # df
        annovar_file = f'{self.outdir}/{self.sample}.{buildver}_multianno.txt'
        df_annovar = utils.parse_annovar(annovar_file)
        return df_annovar


@utils.add_log
def analysis_snp(args):
    step = 'analysis_snp'
    step_snp = Analysis_variant(args, step)
    step_snp.run()


def get_opts_analysis_snp(parser, sub_program):
    parser.add_argument('--annovar_config', help='annovar soft config file.', required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--filter_vcf', help='filter vcf file.', required=True)
        parser.add_argument('--CID_file', help='CID_file.', required=True)
        parser.add_argument('--filter_variant_count_file', help='filter variant count file.', required=True)
