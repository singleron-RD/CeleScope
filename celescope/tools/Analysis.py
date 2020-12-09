import os
import sys
import json
import numpy as np
import pandas as pd
import glob
from celescope.tools.report import reporter
from celescope.tools.utils import glob_genomeDir, log, parse_annovar
from celescope.tools.utils import cluster_tsne_list, parse_match_dir

class Analysis():

    def __init__(
        self,
        sample,
        outdir,
        assay,
        match_dir,
        step,     
    ):
        self.sample = sample
        self.outdir = outdir
        self.assay = assay
        self.match_dir = match_dir
        self.step = step

        match_dict = parse_match_dir(match_dir)
        tsne_df_file = match_dict['tsne_coord']
        marker_df_file = match_dict['markers']
        self.tsne_df = pd.read_csv(tsne_df_file, sep="\t")
        self.marker_df = pd.read_csv(marker_df_file, sep="\t")
        self.tsne_df.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
        self.cluster_tsne = cluster_tsne_list(self.tsne_df)

        if not os.path.exists(outdir):
            os.system('mkdir -p %s' % outdir)

    def add_attrs(self, *kwargs):
        for kwarg in kwargs:
            setattr(self, kwarg, kwargs[kwarg])

    @staticmethod
    def get_table(title, id, df_table):
        """
        return html code
        """
        table_dict = {}
        table_dict['title'] = title
        table_dict['table'] = df_table.to_html(
            escape=False,
            index=False,
            table_id=id,
            justify="center")
        return table_dict

    def report(self, stat=True):
        if stat:
            stat_file = self.outdir + "/stat.txt"
        else:
            stat_file = ''
        t = reporter(
        name=self.step,
        assay=self.assay,
        sample=self.sample,
        outdir=self.outdir + '/..',
        stat_file=stat_file)
        t.get_report()

    def get_cluster_tsne(self, colname, show_tag=True, dfname='tsne_df'):
        """
        tSNE_1	tSNE_2	cluster Gene_Counts
        return data list
        """
        tsne_df = getattr(self, dfname)
        sum_df = tsne_df.groupby([colname]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res = []
        for cluster in sorted(tsne_df[colname].unique()):
            sub_df = tsne_df[tsne_df[colname] == cluster]
            if show_tag:
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

    def report_prepare(self, **kwargs):
        json_file = self.outdir + '/../.data.json'
        if not os.path.exists(json_file):
            data = {}
        else:
            fh = open(json_file)
            data = json.load(fh)
            fh.close()

        for key in kwargs:
            data[key] = kwargs[key]

        with open(json_file, 'w') as fh:
            json.dump(data, fh)

    def get_marker_gene_table(self):
        marker_df = self.process_marker_table()
        table_dict = Analysis.get_table(
            title='Marker Genes by Cluster',
            id='marker_gene_table',
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

    def run(self):
        cluster_tsne = self.get_cluster_tsne(colname='cluster')
        gene_tsne = self.get_gene_tsne()
        table_dict = self.get_marker_gene_table()
        self.report_prepare(
            cluster_tsne=cluster_tsne,
            gene_tsne=gene_tsne,
            table_dict=table_dict,
        )
        self.report()