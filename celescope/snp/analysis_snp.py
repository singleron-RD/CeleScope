
import numpy as np
import pandas as pd
import subprocess
import configparser
import pysam
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from celescope.tools.report import reporter
from celescope.tools.utils import glob_genomeDir, log, parse_annovar
from celescope.tools.utils import parse_vcf
from mutract.utils import read_CID
from celescope.tools.Analysis import Analysis


class analysis_variant(Analysis):
        
    def get_df_count_tsne(self, variant_count_file, CID_file):
        df_vc = pd.read_csv(variant_count_file, sep='\t')
        df_vc = df_vc[df_vc["alt_count"] > 0]
        df_vc_cell = df_vc.groupby('CID').agg({
            'alt_count':'count',
            'VID':list,
        })

        df_CID, df_valid = read_CID(CID_file)
        df_CID = df_CID.reset_index()
        self.tsne_df = pd.merge(self.tsne_df, df_CID, on='barcode', how='left')

        df_vc_barcode = pd.merge(df_vc_cell, df_CID, on='CID')
        df_vc_barcode_tsne = pd.merge(df_vc_barcode, self.tsne_df, on=['barcode','CID'], how='right')
        df_vc_barcode_tsne['value'] = df_vc_barcode_tsne['alt_count']
        df_vc_barcode_tsne['value'] = df_vc_barcode_tsne['value'].fillna(0)
        df_vc_barcode_tsne['value'].astype('int32')
        self.df_count_tsne = df_vc_barcode_tsne

    def get_count_tsne(self):
        def return_text(row):
            text = f'CID:{row["CID"]} <br>Variants:{str(int(row["value"]))} <br>VID:{row["VID"]}'
            return text
        tSNE_1 = list(self.df_count_tsne.tSNE_1)
        tSNE_2 = list(self.df_count_tsne.tSNE_2)
        print(self.df_count_tsne)
        text = list(self.df_count_tsne.apply(return_text, axis=1))
        value = list(self.df_count_tsne.value)
        title = 't-SNE plot Colored by Cell Variant Counts'
        count_tsne = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "text": text, 'value':value, 'title':title}
        return count_tsne

    def add_GT(self, vcf_file):
        vcf = pysam.VariantFile(vcf_file, 'r')
        out_vcf_file = f'{self.outdir}/{self.sample}_addGT.vcf'
        out_vcf = pysam.VariantFile(out_vcf_file, 'w', header=vcf.header)
        for rec in vcf:
            for sample in rec.samples:
                rec.samples[sample]["GT"] = (1,1)
                out_vcf.write(rec)
        vcf.close()
        out_vcf.close()
        return out_vcf_file


    def get_df_table(self, vcf_file, annovar_config):
        
        vcf_file = self.add_GT(vcf_file)

        df_vcf = parse_vcf(vcf_file, infos=['VID','CID'])
        df_annovar = self.annovar(vcf_file, annovar_config)
        df_vcf = pd.concat((df_vcf, df_annovar), axis=1)
        df_vcf["nCell"] = df_vcf["CID"].apply(func=lambda row:1 if isinstance(row,str) else len(row))

        out_df_vcf = f'{self.outdir}/{self.sample}_variant_table.tsv'
        df_vcf.to_csv(out_df_vcf, sep='\t', index=False)

        cols = ['VID','Chrom','Pos','Alleles','Gene','nCell','mRNA','Protein','COSMIC']
        df_vcf = df_vcf[cols]
        return df_vcf


    def run(self, variant_count_file, CID_file, vcf_file, annovar_config):
        cluster_tsne = self.get_cluster_tsne(colname='cluster')

        self.get_df_count_tsne(variant_count_file, CID_file)
        count_tsne = self.get_count_tsne()

        df_vcf = self.get_df_table(vcf_file, annovar_config)
        table_dict = Analysis.get_table(title='Variant table', id='variant_table', df_table=df_vcf)

        self.report_prepare(
            cluster_tsne=cluster_tsne,
            count_tsne=count_tsne,
            table_dict=table_dict,
        )
        self.report(stat=False)

    @log
    def annovar(self, vcf_file, annovar_config):

        # config
        config = configparser.ConfigParser()
        config.read(annovar_config)
        section = config['ANNOVAR']
        dir = section['dir']
        db = section['db']
        buildver = section['buildver']
        protocol = section['protocol']
        operation = section['operation']

        # convert
        input_file = f'{self.outdir}/{self.sample}.input'
        cmd = (
            f'perl {dir}/convert2annovar.pl '
            f'-format vcf4 '
            f'--includeinfo '
            f'{vcf_file} > {input_file}'
        )
        subprocess.check_call(cmd, shell=True)

        # annotate
        cmd = (
            f'perl {dir}/table_annovar.pl '
            f'{input_file} '
            f'{db} '
            f'-buildver {buildver} '
            f'-protocol {protocol} '
            f'-operation {operation} '
            f'-out {self.outdir}/{self.sample} '
            f'--otherinfo '
        )
        analysis_variant.annovar.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

        # df
        annovar_file = f'{self.outdir}/{self.sample}.{buildver}_multianno.txt'
        df_annovar = parse_annovar(annovar_file)
        return df_annovar


@log
def analysis_snp(args):
    step = 'analysis_snp'
    step_snp = analysis_variant(
        args.sample,
        args.outdir,
        args.assay,
        args.match_dir,
        step, 
    )
    step_snp.run(args.variant_count_file, args.CID_file, args.vcf, args.annovar_config)

def get_opts_analysis_snp(parser, sub_program):
    parser.add_argument('--annovar_config', help='annovar soft config file', required=True)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
        parser.add_argument('--vcf', help='vcf file', required=True)
        parser.add_argument('--CID_file', help='CID_file', required=True)
        parser.add_argument('--variant_count_file', help='variant count file', required=True)
        parser.add_argument('--assay', help='assay', required=True)