import subprocess

import pandas as pd

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH
from celescope.rna.mkref import parse_genomeDir_rna
from celescope.tools.step import Step


class AnalysisMixin():
    """
    mixin class for analysis
    child class must inherite Step class
    """

    def __init__(self, args):
        self.args = args
        self.match_dir = None
        if hasattr(args, "match_dir") and args.match_dir:
            self.match_dir = args.match_dir
            self.read_match_dir()
        elif hasattr(args, "tsne_file") and args.tsne_file:
            tsne_df_file = args.tsne_file
            self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
            self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            marker_df_file = args.marker_file
            self.marker_df = pd.read_csv(marker_df_file, sep="\t")
        else:
            self.match_dir = args.outdir + "/../"  # use self

    @utils.add_log
    def seurat(self, matrix_file, save_rds, genomeDir):
        app = ROOT_PATH + "/tools/run_analysis.R"
        genome = parse_genomeDir_rna(genomeDir)
        mt_gene_list = genome['mt_gene_list']
        cmd = (
            f'Rscript {app} '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
            f'--matrix_file {matrix_file} '
            f'--mt_gene_list {mt_gene_list} '
            f'--save_rds {save_rds}'
        )
        self.seurat.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def auto_assign(self, type_marker_tsv):
        rds = f'{self.outdir}/{self.sample}.rds'
        app = ROOT_PATH + "/tools/auto_assign.R"
        cmd = (
            f'Rscript {app} '
            f'--rds {rds} '
            f'--type_marker_tsv {type_marker_tsv} '
            f'--outdir {self.outdir} '
            f'--sample {self.sample} '
        )
        self.auto_assign.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    def get_cluster_tsne(colname, tsne_df, show_colname=True):
        """
        tSNE_1	tSNE_2	cluster Gene_Counts
        return data list
        """

        sum_df = tsne_df.groupby([colname]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res = []
        for cluster in sorted(tsne_df[colname].unique()):
            sub_df = tsne_df[tsne_df[colname] == cluster]
            if show_colname:
                name = f"{colname} {cluster}({percent_df[cluster]}%)"
            else:
                name = f"{cluster}({percent_df[cluster]}%)"
            tSNE_1 = list(sub_df.tSNE_1)
            tSNE_2 = list(sub_df.tSNE_2)
            res.append({"name": name, "tSNE_1": tSNE_1, "tSNE_2": tSNE_2})
        return res

    @staticmethod
    def get_table_dict(marker_table):
        table_dict = Step.get_table(
            title='Marker Genes by Cluster',
            table_id='marker_gene_table',
            df_table=marker_table,
        )
        return table_dict

    def get_marker_table(self):
        """
        return html code
        """

        avg_logfc_col = "avg_log2FC"  # seurat 4
        if "avg_logFC" in self.marker_df.columns:  # seurat 2.3.4
            avg_logfc_col = "avg_logFC"
        marker_df = self.marker_df.loc[:,
                                       ["cluster", "gene", avg_logfc_col, "pct.1", "pct.2", "p_val_adj"]
                                       ]
        marker_df["cluster"] = marker_df["cluster"].apply(lambda x: f"cluster {x}")

        return marker_df

    def get_gene_tsne(self):
        """
        return data dic
        """
        tsne_df = self.tsne_df
        tSNE_1 = list(tsne_df.tSNE_1)
        tSNE_2 = list(tsne_df.tSNE_2)
        Gene_Counts = list(tsne_df.Gene_Counts)
        res = {"tSNE_1": tSNE_1, "tSNE_2": tSNE_2, "Gene_Counts": Gene_Counts}
        return res

    def read_match_dir(self):
        """
        if match_dir is not self, should read match_dir at init
        if it is self, read at run_analysis - need to run seurat first
        """
        if self.match_dir:
            match_dict = utils.parse_match_dir(self.match_dir)
            tsne_df_file = match_dict['tsne_coord']
            self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
            self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            marker_df_file = match_dict['markers']
            self.marker_df = pd.read_csv(marker_df_file, sep="\t")

    def run_analysis(self):
        self.read_match_dir()
        self.cluster_tsne = self.get_cluster_tsne(colname='cluster', tsne_df=self.tsne_df)
        self.gene_tsne = self.get_gene_tsne()
        marker_table = self.get_marker_table()
        self.table_dict = self.get_table_dict(marker_table)
