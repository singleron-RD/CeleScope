import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT, ROOT_PATH


def get_opts_analysis_cite(parser, sub_program):
    if sub_program:
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--citeseq_mtx', help='citeseq matrix .gz file', required=True)
        s_common(parser)


def analysis_cite(args):
    with Analysis_cite(args, display_title='Analysis') as runner:
        runner.run()


class Analysis_cite(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        # data
        self.tsne_dict = utils.parse_match_dir(args.match_dir)

    def run(self):

        rds = self.tsne_dict['rds']
        app = ROOT_PATH + "/citeseq/analysis_cite.R"
        cmd = (
            f'Rscript {app} '
            f'--rds {rds} '
            f'--citeseq_mtx {self.args.citeseq_mtx} '
            f'--outdir {self.args.outdir} '
            f'--sample {self.args.sample} '
        )
        self.debug_subprocess_call(cmd)
