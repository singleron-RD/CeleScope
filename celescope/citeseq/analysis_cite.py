import subprocess

from celescope.tools import utils
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
        self.rds = f"{self.out_prefix}.rds"
        match_dict = utils.parse_match_dir(args.match_dir)
        self.match_matrix_dir = match_dict['matrix_dir']

    def run_citeseq(self):

        app = ROOT_PATH + "/citeseq/analysis_cite.R"
        cmd = (
            f'Rscript {app} '
            f'--rds {self.rds} '
            f'--citeseq_mtx {self.args.citeseq_mtx} '
            f'--outdir {self.outdir} '
            f'--sample {self.sample} '
            '2>&1 '
        )
        self.debug_subprocess_call(cmd)
    
    @utils.add_log
    def run_seurat(self):
        app = ROOT_PATH + "/citeseq/run_analysis.R"
        cmd = (
            f'Rscript {app} '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
            f'--matrix_file {self.match_matrix_dir} '
            f'--mt_gene_list None '
            f'--save_rds True '
            '2>&1 '
        )
        self.run_seurat.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run(self):
        self.run_seurat()
        self.run_citeseq()
