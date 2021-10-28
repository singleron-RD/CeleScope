import configparser
import subprocess

import pandas as pd
import pysam
from venn import generate_petal_labels,draw_venn,generate_colors

import celescope.tools.utils as utils
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common
from celescope.snp.variant_calling import read_CID
from celescope.__init__ import HELP_DICT


class Analysis_variant(Step, AnalysisMixin):
    """
    Features
    - Annotate variants with [Annovar](https://annovar.openbioinformatics.org/en/latest/).

    Output
    `{sample}.{genome_version}_multianno.txt` Annovar main output file. `CID` and `VID` are added to the `Otherinfo` column.

    `{sample}_variant_table.tsv` Formatted `multianno` file with `ncell_cover`and `ncell_alt` added.
    
    `{sample}_variant_top5.jpg` The Venn diagram of the 5 variants with the highest `ncell_alt`.

    `{sample}_variant_ncell.tsv` Number of cells with read count at each variant's position. 
    - `VID`: Variant ID. 
    - `ncell_cover`: number of cells with read count at this position. 
    - `ncell_alt`: number of cells with variant read count only. 
    - `ncell_ref`: number of cells with reference read count only. 
    - `ncell_ref_and_alt`: number of cells with both variant and reference read count.
    - `RID`: Target region ID. This column will be added when `--panel` option were provided.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        AnalysisMixin.__init__(self, args)
        self.filter_variant_count_file = args.filter_variant_count_file
        self.CID_file = args.CID_file
        self.filter_vcf = args.filter_vcf
        self.annovar_config = args.annovar_config
        self.match_dir = args.match_dir

        # parse
        self.annovar_section = self.read_annovar_config()

        # out
        self.ncell_file = f'{self.out_prefix}_variant_ncell.tsv'
        self.vcf_GT = f'{self.out_prefix}_addGT.vcf'
        buildver = self.annovar_section['buildver']
        self.annovar_file = f'{self.out_prefix}.{buildver}_multianno.txt'

    @utils.add_log
    def ncell_metrics(self):
        variant_count_df = pd.read_table(self.filter_variant_count_file, sep='\t')
        variant_count_df = variant_count_df.loc[(variant_count_df["ref_count"]!=0) | 
            (variant_count_df["alt_count"]!=0)]
        if 'RID' in variant_count_df.columns: 
            df_cover = variant_count_df.groupby('VID').agg({'RID':'first','CID':'count'})
        else:
            df_cover = variant_count_df.groupby('VID').agg({'CID': 'count'})
        df_cover.rename(columns={'CID': 'ncell_cover'}, inplace=True)

        df_ref = variant_count_df.loc[(variant_count_df["ref_count"]!=0) & 
            (variant_count_df["alt_count"]==0)].groupby('VID').agg({'CID':'count'})
        df_ref.rename(columns={'CID': 'ncell_ref'}, inplace=True)

        df_alt = variant_count_df.loc[(variant_count_df["ref_count"]==0) & 
            (variant_count_df["alt_count"]!=0)].groupby('VID').agg({'CID':'count'})
        df_alt.rename(columns={'CID': 'ncell_alt'}, inplace=True)

        df_ref_and_alt = variant_count_df.loc[(variant_count_df["ref_count"]!=0) & 
            (variant_count_df["alt_count"]!=0)].groupby('VID').agg({'CID':'count'})
        df_ref_and_alt.rename(columns={'CID': 'ncell_ref_and_alt'}, inplace=True)

        df_list = [df_ref, df_alt, df_ref_and_alt]
        df_ncell = df_cover
        for df in df_list:
            df_ncell = pd.merge(df_ncell, df, on='VID', how='left')
        df_ncell.fillna(0, inplace=True)

        # dtype after fillna is float
        float_cols = ['ncell_ref', 'ncell_alt', 'ncell_ref_and_alt']
        for col in float_cols:
            df_ncell[col] = df_ncell[col].astype('int')

        df_ncell.to_csv(self.ncell_file, sep = '\t',index = True)

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
        with pysam.VariantFile(self.filter_vcf, 'r') as vcf:
            with pysam.VariantFile(self.vcf_GT, 'w', header=vcf.header) as vcf_GT:
                for rec in vcf:
                    for sample in rec.samples:
                        rec.samples[sample]["GT"] = (1, 1)
                        vcf_GT.write(rec)

    def get_df_table(self):

        df_vcf = utils.parse_vcf(self.vcf_GT, infos=['VID', 'CID'])
        df_annovar = self.annovar()
        df_vcf = pd.concat((df_vcf, df_annovar), axis=1)
        ncell_df = pd.read_table(self.ncell_file,sep = "\t")
        ncell_df.loc[:,"VID"] = ncell_df.loc[:,"VID"].astype(str)
        df_vcf["nCell"] = df_vcf["CID"].apply(func=lambda row: 1 if isinstance(row, str) else len(row))

        df_vcf = pd.merge(left = df_vcf,
                          right = ncell_df,
                          on = "VID",
                          how = "left")

        out_df_vcf = f'{self.outdir}/{self.sample}_variant_table.tsv'
        df_vcf.drop("nCell",axis = 1).to_csv(out_df_vcf, sep='\t', index=False)

        cols = ['VID', "CID",'Chrom', 'Pos', 'Alleles', 'Gene',  
            'ncell_cover', "ncell_ref", 'ncell_alt', "ncell_ref_and_alt",
            'mRNA', 'Protein', 'COSMIC']
        df_vcf = df_vcf[cols]
        return df_vcf
    
    def get_venn_plot(self):
        df_top_5 = self.get_df_table().sort_values(by = "ncell_alt",ascending=False).iloc[:5,:]
        plot = {}
        cid_lst = df_top_5.loc[:,"CID"].to_list()
        vid_lst = df_top_5.loc[:,"VID"].to_list()
        for cid,vid in zip(cid_lst,vid_lst):
            plot[f"VID_{vid}"] = set(cid)
        share_cid  = list(set.intersection(*map(set,cid_lst)))
        if share_cid == []:
            share_cid.append("None")
        #venn plot
        set_cid = list(plot.values())
        set_name = list(plot.keys())
        labels = generate_petal_labels(set_cid)
        plot = draw_venn(
                         petal_labels=labels, 
                         dataset_labels=set_name,
                         hint_hidden=False,
                         colors=generate_colors(n_colors=5), 
                         figsize=(8, 8),
                         fontsize=14, 
                         legend_loc="best",
                         ax=None
                         )
        fig = plot.get_figure()
        fig.savefig(f'{self.outdir}/{self.sample}_variant_top5.jpg',dpi = 600)
        pd.DataFrame({"top5_variant_shared_cells":share_cid}).to_csv(f'{self.outdir}/{self.sample}_top5_shared_cells.tsv',sep = '\t',index = None)

    def run(self):
        self.add_GT()
        self.ncell_metrics()
        cluster_tsne = self.get_cluster_tsne(colname='cluster', tsne_df=self.tsne_df)
        df_count_tsne = self.get_df_count_tsne()
        count_tsne = self.get_count_tsne(df_count_tsne)
        df_vcf = self.get_df_table()
        table_dict = Step.get_table(title='Variant table', table_id='variant_table', df_table=df_vcf.drop(["CID"],axis = 1))

        self.add_data_item(cluster_tsne=cluster_tsne)
        self.add_data_item(count_tsne=count_tsne)
        self.add_data_item(table_dict=table_dict)
        self.clean_up()
        self.get_venn_plot()

    def read_annovar_config(self):
        '''
        read annovar config file
        '''
        config = configparser.ConfigParser()
        config.read(self.annovar_config)
        section = config['ANNOVAR']
        return section

    @utils.add_log
    def annovar(self):

        section = self.annovar_section
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
        df_annovar = utils.parse_annovar(self.annovar_file)
        return df_annovar


@utils.add_log
def analysis_snp(args):
    step = 'analysis_snp'
    step_snp = Analysis_variant(args, step)
    step_snp.run()


def get_opts_analysis_snp(parser, sub_program):
    parser.add_argument('--annovar_config', help='ANNOVAR config file.', required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--filter_vcf', help='filter vcf file.', required=True)
        parser.add_argument('--CID_file', help='CID_file.', required=True)
        parser.add_argument('--filter_variant_count_file', help='filter variant count file.', required=True)

