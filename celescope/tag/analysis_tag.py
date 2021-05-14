import pandas as pd

import celescope.tools.utils as utils
from celescope.tools.analysis_mixin import AnalysisMixin
from celescope.tools.step import Step, s_common


class Analysis_tag(Step, AnalysisMixin):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        AnalysisMixin.__init__(self, args)

    def run(self):
        self.run_analysis()
        tsne_tag_df = pd.read_csv(self.args.tsne_tag_file, sep="\t", index_col=0)
        feature_tsne = self.get_cluster_tsne(colname='tag', tsne_df=tsne_tag_df, show_colname=False)
        self.add_data_item(cluster_tsne=self.cluster_tsne)
        self.add_data_item(feature_tsne=feature_tsne)
        self.add_data_item(table_dict=self.table_dict)

        self.clean_up()


def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--tsne_tag_file', help='tsne tag file', required=True)
        parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path", required=True)


@utils.add_log
def analysis_tag(args):
    step_name = 'analysis_tag'
    ana = Analysis_tag(args, step_name)
    ana.run()