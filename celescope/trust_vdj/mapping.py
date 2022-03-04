import pandas as pd
import numpy as np 
from celescope.tools import utils
from celescope.tools.step import Step, s_common
import celescope
import os
import subprocess
import glob
from celescope.trust_vdj.__init__ import TOOLS_DIR

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
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.seqtype = args.seqtype
        self.match_dir = args.match_dir

        try:
            self.rds = glob.glob(f'{self.match_dir}/06.analysis/*.rds')[0]
            self.assign_file = glob.glob(f'{self.match_dir}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
        except:
            pass
  
        self.contig = glob.glob(f'{self.outdir}/../04.summarize/{self.sample}_filtered_contig.csv')[0]
        self.summary = []

        if self.seqtype == 'TCR':
            self.Celltype = {'T_cells','NKT_cells','T cells','NK T cells','Tcells'}
            self._name = "Tcells"
        elif self.seqtype == 'BCR':
            self.Celltype = {'Plasma_cells','B_cells','Mature_B_cell','Plasma cells','B cells','Bcells'}
            self._name = "Bcells"

    @utils.add_log
    def process(self):
        run_mapping(self.rds,self.contig,self.sample,self.outdir,self.assign_file)
        meta = pd.read_csv(glob.glob(f'{self.outdir}/{self.sample}_meta.csv')[0])
        metaTB = meta[meta['CellTypes'].isin(self.Celltype)]
        mappedmeta = meta[meta['Class']=='T/BCR']
        mappedmetaTB = mappedmeta[mappedmeta['CellTypes'].isin(self.Celltype)]
        
        Transcriptome_cell_number = meta.shape[0]
        TB_cell_number = metaTB.shape[0]
        Mapped_Transcriptome_cell_number = mappedmeta.shape[0]
        Mapped_TB_cell_number = mappedmetaTB.shape[0]

        self.summary.append({
            'item': 'Cell Number in Matched transcriptome',
            'count': Transcriptome_cell_number,
            'total_count': np.nan
        })
        self.summary.append({
            'item': 'Cell Number Successfully Mapped to transcriptome',
            'count': Mapped_Transcriptome_cell_number,
            'total_count': Transcriptome_cell_number
        })
        self.summary.append({
            'item': f'{self._name} Number in Matched transcriptome',
            'count': TB_cell_number,
            'total_count': np.nan
        })
        self.summary.append({
            'item': f'Cell Number Successfully Mapped to {self._name} in transcriptome',
            'count': Mapped_TB_cell_number,
            'total_count': TB_cell_number
        })

        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(self.summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file) 
    
    def run(self):
        try:
            assert self.rds and self.assign_file
            self.process()

        except AssertionError and AttributeError:
            print("rds file and type file do not exist" + "\n" )
        
        except ZeroDivisionError:
            print(f"Not found auto-assigned {self._name} in matched sc-RNA")


def mapping(args):
    step_name = 'mapping'
    mapping_obj = Mapping(args, step_name)
    mapping_obj.run()

    
def get_opts_mapping(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR',
                        choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
    return parser