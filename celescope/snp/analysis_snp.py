import configparser
import subprocess

import pandas as pd
from celescope.capture_virus.analysis_capture_virus import Analysis_capture_virus
import pysam
from venn import generate_petal_labels,draw_venn,generate_colors

import celescope.tools.utils as utils
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT, ROOT_PATH


class Analysis_variant(Step, AnalysisMixin):
    """
    Features
    - Annotate variants with [Annovar](https://annovar.openbioinformatics.org/en/latest/).

    Output


    """

    def __init__(self, args,display_title=None):
        Step.__init__(self, args,display_title=display_title)
        AnalysisMixin.__init__(self, args)
        self.vcf_file = args.vcf
        self.annovar_config = args.annovar_config
        self.match_dir = args.match_dir

        # parse
        self.annovar_section = self.read_annovar_config()
        buildver = self.annovar_section['buildver']

        # out
        self.annovar_outdir = f'{self.outdir}/annovar/'
        utils.check_mkdir(self.annovar_outdir)
        self.input_file = f'{self.annovar_outdir}/{self.sample}.input'
        self.multianno_file = f'{self.annovar_outdir}/{self.sample}.{buildver}_multianno.txt'

        self.gt_file = f'{self.out_prefix}_gt.csv'
        self.ncell_file = f'{self.out_prefix}_variant_ncell.csv'
        self.variant_table_file = f'{self.out_prefix}_variant_table.csv'

    @utils.add_log
    def write_gt(self):
        app = f'{ROOT_PATH}/snp/vcfR.R'
        cmd = (
            f'Rscript {app} '
            f'--vcf {self.vcf_file} '
            f'--out {self.gt_file} '
        )
        self.debug_subprocess_call(cmd)

    @utils.add_log
    def write_ncell(self):
        """
        parse gt_file to collect each genotype cell count into ncell_file
        """
        df = pd.read_csv(self.gt_file, index_col=0)
        df_ncell = df.apply(pd.Series.value_counts, axis=1).fillna(0).astype(int)
        df_ncell.to_csv(self.ncell_file, index = True)

    @utils.add_log
    def run_annovar(self):

        section = self.annovar_section
        annovar_dir = section['dir']
        db = section['db']
        buildver = section['buildver']
        protocol = section['protocol']
        operation = section['operation']

        # convert
        cmd = (
            f'perl {annovar_dir}/convert2annovar.pl '
            f'--format vcf4old '
            f'--includeinfo '
            f'{self.vcf_file} > {self.input_file}'
        )
        self.debug_subprocess_call(cmd)

        # annotate
        cmd = (
            f'perl {annovar_dir}/table_annovar.pl '
            f'{self.input_file} '
            f'{db} '
            f'-buildver {buildver} '
            f'-protocol {protocol} '
            f'-operation {operation} '
            f'-out {self.annovar_outdir}/{self.sample} '
            f'--otherinfo '
        )
        self.debug_subprocess_call(cmd)

    def get_variant_table(self):

        df_vcf = utils.parse_vcf(self.vcf_file, infos=[])
        df_annovar = utils.parse_annovar(self.multianno_file)
        df_vcf = pd.concat((df_vcf, df_annovar), axis=1)
        df_ncell = pd.read_csv(self.ncell_file)
        df_vcf = pd.concat([df_vcf, df_ncell], axis=1)

        cols = ['Chrom', 'Pos', 'Alleles', 'Gene', '0/0', "0/1", '1/1', 'mRNA', 'Protein', 'COSMIC']
        df_vcf = df_vcf[cols]
        self.variant_table = df_vcf
        self.variant_table.to_csv(self.variant_table_file, index=False)
    
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
        self.write_gt()
        self.write_ncell()
        self.run_annovar()
        self.get_variant_table()
        table_dict = self.get_table(title='Variant table', table_id='variant_table', df_table=self.variant_table)
        self.add_data_item(table_dict=table_dict)
        self.clean_up()
        #self.get_venn_plot()

    def read_annovar_config(self):
        '''
        read annovar config file
        '''
        config = configparser.ConfigParser()
        config.read(self.annovar_config)
        section = config['ANNOVAR']
        return section


@utils.add_log
def analysis_snp(args):
    with Analysis_variant(args,display_title="Analysis") as runner:
        runner.run()


def get_opts_analysis_snp(parser, sub_program):
    parser.add_argument('--annovar_config', help='ANNOVAR config file.', required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--vcf', help='vcf file.', required=True)

