import numpy as np

from celescope.tools import utils, analysis_wrapper, plotly_plot
from celescope.tools.tag.analysis_tag import Analysis_tag as At
from celescope.tools.step import s_common
from celescope.__init__ import HELP_DICT

# default colnames in tsne_tag file 
DEFAULT_COLS = ['', 'tSNE_1', 'tSNE_2', 'cluster', 'Gene_Counts']

class Analysis_tag(At):
    '''
    ## Features
    - Combine scRNA-Seq clustering infromation with tag assignment.
    '''
    def run(self):
        colnames = list(self.df_tsne_tag.columns)
        tag_cols = set(colnames) - set(DEFAULT_COLS)
        if len(tag_cols) > 1:
            raise ValueError(f"Error: More than one tag barcode in tag_barcode_fasta: {tag_cols}. Currently, Multiple tag barcodes are not supported to display in HTML report.")
        tag_col = list(tag_cols)[0]

        log_tag_col = tag_col + '_log2'
        self.df_tsne_tag[log_tag_col] = np.log2(self.df_tsne_tag[tag_col] + 1)

        report_runner = analysis_wrapper.Report_runner(self.args, display_title=self.display_title)
        df_tsne, _df_marker = report_runner.get_df()

        tsne_cluster = plotly_plot.Tsne_plot(df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_tag = plotly_plot.Tsne_plot(self.df_tsne_tag, log_tag_col, discrete=False).get_plotly_div()
        self.add_data(tsne_tag=tsne_tag)

def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser.add_argument('--tsne_tag_file', help='`{sample}_tsne_tag.tsv` from count_tag. ', required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--tsne_file", help=HELP_DICT['tsne_file'])
        parser = s_common(parser)

@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args, display_title="Analysis") as runner:
        runner.run()

