import pandas as pd
import pysam
import subprocess
import os
from venn import generate_petal_labels, draw_venn, generate_colors

from celescope.tools import utils
from celescope.tools.step import Step
from celescope.tools.step import s_common
from celescope.tools.target_metrics import get_gene_list
from celescope.__init__ import HELP_DICT, ROOT_PATH
from celescope.snp.__init__ import PANEL


AA_DICT = {
    'Gly' : 'G',
    'Ala' : 'A',
    'Val' : 'V',
    'Leu' : 'L',
    'Ile' : 'I',
    'Phe' : 'F',
    'Trp' : 'W',
    'Tyr' : 'Y',
    'Asp' : 'D',
    'Asn' : 'N',
    'Glu' : 'E',
    'Lys' : 'K',
    'Gln' : 'Q',
    'Met' : 'M',
    'Ser' : 'S',
    'Thr' : 'T',
    'Cys' : 'C',
    'Pro' : 'P',
    'His' : 'H',
    'Arg' : 'R',
}


def parse_variant_ann(variant_ann_file):
    """
    Args:
        variant_ann_file: variant annotation file from snpEff.
    """
    gene_list, mRNA_list, protein_list = [], [], []

    with open(variant_ann_file) as f:
        for line in f.readlines():
            if not line.startswith("#"):
                info = line.split('\t')[7]
                anns = info.split("|")
                gene = anns[3]
                gene_list.append(gene)
            
                tmp1, tmp2 = [], []
                for ann in anns:
                    if ann.startswith("c."):
                        exon_loc = anns[anns.index(ann) - 1].split('/')[0]
                        exon = ann.strip("c.")
                        exon = f"exon{exon_loc}:{exon}"
                        if exon not in tmp1:
                            tmp1.append(exon)

                    if ann.startswith("p."):
                        protein = ann.strip("p.")
                        for i in AA_DICT:
                            protein = protein.replace(i, AA_DICT[i])
                        if protein not in tmp2:
                            tmp2.append(protein)
                        
                mRNA_list.append(','.join(tmp1))
                protein_list.append(','.join(tmp2))

    return (gene_list, mRNA_list, protein_list)


def parse_vcf_to_df(vcf_file, cols=('chrom', 'pos', 'alleles'), infos=('VID', 'CID')):
    """
    Read cols and infos into pandas df
    """
    vcf = pysam.VariantFile(vcf_file)
    df = pd.DataFrame(columns=[col.capitalize() for col in cols] + infos)
    rec_dict = {}
    for rec in vcf.fetch():

        for col in cols:
            rec_dict[col.capitalize()] = getattr(rec, col)
            if col == 'alleles':
                rec_dict['Alleles'] = '-'.join(rec_dict['Alleles'])

        for info in infos:
            rec_dict[info] = rec.info[info]

        '''
        rec_dict['GT'] = [s['GT'] for s in rec.samples.values()][0]
        rec_dict['GT'] = [str(item) for item in rec_dict['GT']]
        rec_dict['GT'] = '/'.join(rec_dict['GT'])
        '''
        df_new = pd.DataFrame(rec_dict, index=[0])
        df = pd.concat([df, df_new])

    vcf.close()
    df.reset_index(drop=True, inplace=True)
    return df


class Analysis_snp(Step):
    """
    ## Features
    - Annotate variants with [snpEff](http://pcingola.github.io/SnpEff/).

    ## Output
    - `{sample}_gt.csv` Genotypes of variants of each cell. Rows are variants and columns are cells.
    - `{sample}_variant_ncell.csv` Number of cells with each genotype.
    - `{sample}_variant_table.csv` annotated with snpEff.

    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.vcf_file = args.vcf

        # parse
        self.gene_list, self.n_gene = get_gene_list(args)

        # data
        self.variant_table = None

        # out
        self.snpeff_outdir = f'{self.outdir}/snpEff/'
        self.snpeff_ann = f'{self.snpeff_outdir}/variants_ann.vcf'
        utils.check_mkdir(self.snpeff_outdir)

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
            '2>&1 '
        )
        self.debug_subprocess_call(cmd)

    @utils.add_log
    def write_ncell(self):
        """
        parse gt_file to collect each genotype cell count into ncell_file
        """
        df = pd.read_csv(self.gt_file, index_col=0)
        df_ncell = df.apply(pd.Series.value_counts, axis=1).fillna(0).astype(int)
        df_ncell.to_csv(self.ncell_file, index=True)

    @utils.add_log
    def run_snpEff(self):
        # Filter -no-downstream -no-upstream -no-utr -no-intron -no-intergenic -no SPLICE_SITE_REGION
        cmd = (
            f"snpEff -Xmx8g -v GRCh38.99 {os.path.abspath(self.vcf_file)} > variants_ann.vcf "
            "-no-downstream -no-upstream -no-utr -no-intron -no-intergenic -no SPLICE_SITE_REGION "
        )
        self.run_snpEff.logger.info(cmd)

        cwd = os.getcwd()
        os.chdir(self.snpeff_outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid can not find '09.analysis_snp/stat.txt' error
        os.chdir(cwd)

    def get_variant_table(self):

        df_vcf = parse_vcf_to_df(self.vcf_file, infos=[])
        ann_result = parse_variant_ann(self.snpeff_ann)
        df_vcf["Gene"], df_vcf["mRNA"], df_vcf["Protein"] = ann_result[0], ann_result[1], ann_result[2]
        df_ncell = pd.read_csv(self.ncell_file)
        df_vcf = pd.concat([df_vcf, df_ncell], axis=1)

        cols = ["Chrom", "Pos", "Alleles", "Gene", "0/0", "0/1", "1/1", "mRNA", "Protein"]
        df_vcf = df_vcf[cols]
        df_vcf = df_vcf[df_vcf.Gene.isin(self.gene_list)]

        self.variant_table = df_vcf
        self.variant_table.reset_index(drop=True, inplace=True)
        self.variant_table.to_csv(self.variant_table_file, index=False)
    

    def get_venn_plot(self):
        df_top_5 = self.get_df_table().sort_values(by="ncell_alt", ascending=False).iloc[:5, :]
        plot = {}
        cid_lst = df_top_5.loc[:, "CID"].to_list()
        vid_lst = df_top_5.loc[:, "VID"].to_list()
        for cid, vid in zip(cid_lst, vid_lst):
            plot[f"VID_{vid}"] = set(cid)
        share_cid = list(set.intersection(*map(set, cid_lst)))
        if share_cid == []:
            share_cid.append("None")
        # venn plot
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
        fig.savefig(f'{self.outdir}/{self.sample}_variant_top5.jpg', dpi=600)
        pd.DataFrame({"top5_variant_shared_cells": share_cid}).to_csv(
            f'{self.outdir}/{self.sample}_top5_shared_cells.tsv', sep='\t', index=None)

    def add_help(self):
        '''
            <p> Chrom : chromosome name.</p>
            <p> Pos : the 1-based position of the variation on the given sequence..</p>
            <p> Alleles : REF(reference base or bases in the case of an indel) - ALT(alternative alleles).</p>
            <p> 0/0, 0/1, 1/1: number of cells with each genotype.</p>
            <p> Gene : gene symbol.</p>
            <p> mRNA :  A standard nomenclature is used in specifying the sequence changes.</p>
            <p> Protein :  A standard nomenclature is used in specifying the sequence changes.</p>
        '''
        self.add_help_content(
            name='Chrom',
            content='Chromosome name'
        )
        self.add_help_content(
            name='Pos',
            content='the 1-based position of the variation on the given sequence'
        )
        self.add_help_content(
            name='Alleles',
            content='REF(reference base or bases in the case of an indel) - ALT(alternative alleles)'
        )
        self.add_help_content(
            name='0/0, 0/1, 1/1',
            content='number of cells with each genotype'
        )
        self.add_help_content(
            name='Gene',
            content='gene symbol'
        )
        self.add_help_content(
            name='mRNA',
            content='A standard nomenclature is used in specifying the sequence changes'
        )
        self.add_help_content(
            name='Protein',
            content='A standard nomenclature is used in specifying the sequence changes'
        )

    def run(self):
        self.write_gt()
        self.write_ncell()
        self.run_snpEff()
        self.get_variant_table()
        self.add_help()
        table_dict = self.get_table_dict(title='Variant table', table_id='variant', df_table=self.variant_table)
        self.add_data(table_dict=table_dict)
        # self.get_venn_plot()


@utils.add_log
def analysis_snp(args):
    with Analysis_snp(args, display_title='Analysis') as runner:
        runner.run()


def get_opts_analysis_snp(parser, sub_program):
    parser.add_argument("--gene_list", help=HELP_DICT['gene_list'])
    parser.add_argument("--panel", help=HELP_DICT['panel'], choices=list(PANEL))
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--vcf', help='vcf file.', required=True)
