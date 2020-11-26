import os
import sys
import json
import numpy as np
import pandas as pd
import glob
from celescope.tools.report import reporter
from celescope.tools.utils import log
from celescope.tools.Analysis import Analysis

class Analysis_tag(Analysis):

    def run(self, tsne_tag_file):
        cluster_tsne = self.get_cluster_tsne(colname='cluster')
        self.tsne_tag_df = pd.read_csv(tsne_tag_file, sep="\t", index_col=0)
        feature_tsne = self.get_cluster_tsne(colname='tag', show_tag=False, dfname='tsne_tag_df')
        marker_gene_table = self.marker_table()
        self.report_prepare(
            cluster_tsne=cluster_tsne, 
            feature_tsne=feature_tsne,
            marker_gene_table=marker_gene_table,
        )
        self.report(stat=False)





