import pandas as pd
from celescope.tools import utils
from celescope.tools.step import Step, s_common
import celescope
import os
import subprocess
import glob

TOOLS_DIR = os.path.dirname(celescope.tools.__file__)

def run_mapping(rds, contig, sample, outdir, assign):
    cmd = (
        f'Rscript {TOOLS_DIR}/VDJmapping.R '
        f'--rds {rds} '
        f'--VDJ {contig} '
        f'--sample {sample} '
        f'--outdir {outdir} '
        f'--assign_file {assign}'
    )
    subprocess.check_call(cmd, shell=True)

class Mapping(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.match_dir = args.match_dir
        self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')[0]
        self.assign_file =  glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
        self.contig = glob.glob(f'{self.outdir}/../05.match/match_contigs.csv')[0]

    @utils.add_log
    def process(self):
        run_mapping(self.rds,self.contig,self.sample,self.outdir,self.assign_file)
    
    def run(self):    
        self.process()

def mapping(args):

    with Mapping(args, display_title="Mapping") as runner:
        runner.run()

    
def get_opts_mapping(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
    return parser