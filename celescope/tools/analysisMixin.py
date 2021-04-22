import os
import sys
import json
import numpy as np
import pandas as pd
import celescope.tools.utils as utils
from celescope.tools.Step import Step


class AnalysisMixin():

    def __init__(self, args):
        if hasattr(args, "match_dir") and args.match_dir:
            self.match_dir = args.match_dir
        else:
            self.match_dir = args.outdir + "/../" # use self

        if self.match_dir:
            match_dict = utils.parse_match_dir(self.match_dir)
            tsne_df_file = match_dict['tsne_coord']
            marker_df_file = match_dict['markers']
            self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
            self.marker_df = pd.read_csv(marker_df_file, sep="\t")
            self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            self.cluster_tsne = utils.cluster_tsne_list(self.tsne_df)

    @staticmethod
    def cluster_tsne_list(tsne_df):
        """
        tSNE_1	tSNE_2	cluster Gene_Counts
        return data list
        """
        sum_df = tsne_df.groupby(["cluster"]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res = []
        for cluster in sorted(tsne_df.cluster.unique()):
            sub_df = tsne_df[tsne_df.cluster == cluster]
            name = "cluster {cluster}({percent}%)".format(
                cluster=cluster, percent=percent_df[cluster])
            tSNE_1 = list(sub_df.tSNE_1)
            tSNE_2 = list(sub_df.tSNE_2)
            res.append({"name": name, "tSNE_1": tSNE_1, "tSNE_2": tSNE_2})
        return res

    def get_cluster_tsne(self, colname, tsne_df, show_colname=True):
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

    def process_marker_table(self):
        """
        return html code
        """
        marker_df = self.marker_df.loc[:, ["cluster", "gene",
                                    "avg_logFC", "pct.1", "pct.2", "p_val_adj"]]
        marker_df["cluster"] = marker_df["cluster"].apply(lambda x: f"cluster {x}")
        return marker_df

    def get_marker_gene_table(self):
        marker_df = self.process_marker_table()
        table_dict = Step.get_table(
            title='Marker Genes by Cluster',
            table_id='marker_gene_table',
            df_table=marker_df,
        )
        return table_dict

    def marker_table(self):
        """
        return html code
        """
        marker_df = self.marker_df.loc[:, ["cluster", "gene",
                                    "avg_logFC", "pct.1", "pct.2", "p_val_adj"]]
        marker_df["cluster"] = marker_df["cluster"].apply(lambda x: f"cluster {x}")
        marker_gene_table = marker_df.to_html(
            escape=False,
            index=False,
            table_id="marker_gene_table",
            justify="center")
        return marker_gene_table

    def run_analysis(self):
        self.cluster_tsne = self.get_cluster_tsne(colname='cluster', tsne_df=self.tsne_df)
        self.gene_tsne = self.get_gene_tsne()
        self.table_dict = self.get_marker_gene_table()