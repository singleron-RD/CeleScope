import pandas as pd
from celescope.snp.analysis_snp import Analysis_variant

import celescope.tools.utils as utils
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common


class Analysis_tag(Step, AnalysisMixin):
    """
    Features
    - Combine scRNA-Seq clustering infromation with tag assignment.
    """

    def __init__(self, args,display_title=None):
        Step.__init__(self, args,display_title=display_title)
        AnalysisMixin.__init__(self, args)

    def add_help(self):
        self.add_help_content(
            name='Marker Genes by Cluster',
            content='differential expression analysis based on the non-parameteric Wilcoxon rank sum test'
        )
        self.add_help_content(
            name='avg_log2FC',
            content='log fold-change of the average expression between the cluster and the rest of the sample'
        )
        self.add_help_content(
            name='pct.1',
            content='The percentage of cells where the gene is detected in the cluster'
        )
        self.add_help_content(
            name='pct.2',
            content='The percentage of cells where the gene is detected in the rest of the sample'
        )
        self.add_help_content(
            name='p_val_adj',
            content='Adjusted p-value, based on bonferroni correction using all genes in the dataset'
        )

    def run(self):
        self.run_analysis()
        tsne_tag_df = pd.read_csv(self.args.tsne_tag_file, sep="\t", index_col=0)
        feature_tsne = self.get_cluster_tsne(colname='tag', tsne_df=tsne_tag_df, show_colname=False)
        self.add_data(cluster_tsne=self.cluster_tsne)
        self.add_data(feature_tsne=feature_tsne)
        self.add_data(table_dict=self.table_dict)
        self.add_help()


def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser.add_argument('--tsne_tag_file', help='`{sample}_tsne_tag.tsv` from count_tag. ', required=True)
        parser.add_argument("--match_dir", help="Match celescope scRNA-Seq directory. ")
        parser.add_argument("--tsne_file", help="t-SNE coord file.")
        parser.add_argument("--marker_file", help="Marker file.")
        parser = s_common(parser)


@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args,display_title="Analysis") as runner:
        runner.run()
