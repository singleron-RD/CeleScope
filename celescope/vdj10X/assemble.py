import os
import subprocess

import numpy as np
import pandas as pd
from celescope.tools import utils
from celescope.tools.cellranger3 import get_plot_elements
from celescope.tools.step import Step, s_common
from celescope.vdj10X.__init__ import ref_dict, soft_dict

pd.set_option('max_colwidth',200)

class Assemble(Step):
    
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        
        self.mem = args.mem
        self.species = args.species
        self.soft = args.soft
        self.thread = args.thread
        self.fqs_dir = args.fqs_dir
        self.count_file = args.count_file
        self.barcode_dic = args.barcode_dic
        
        self.cwd_path = os.getcwd()
        self.outs = f'{self.outdir}/{self.sample}/outs'
        
        
    @utils.add_log
    def run_assemble(self):

        ref_path = ref_dict[self.soft][self.species]
        soft_path = soft_dict[self.soft]
        #cwd_path = os.getcwd()
        cmd = (
            f'{soft_path} vdj '
            f'--id={self.sample} '
            f'--reference={ref_path} '
            f'--fastqs={self.cwd_path}/{self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )
        Assemble.run_assemble.logger.info(cmd)
        with open(f'{self.outdir}/{self.sample}_vdj_10X.sh', 'w') as f:
            f.write(cmd)
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
    
    @utils.add_log   
    def gen_stat_file(self):
        os.chdir(self.cwd_path)
        
        # gen clone table
        clonetype_count = pd.read_csv(f'{self.outs}/clonotypes.csv', sep=',', index_col=None)
        
        clonetype_count['ClonotypeID'] = clonetype_count['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        clonetype_count['Frequency'] = clonetype_count['frequency']
        clonetype_count['Proportion'] = clonetype_count['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        clonetype_count['CDR3_aa'] = clonetype_count['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))
        
        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', clonetype_count[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']])
        self.add_data_item(table_dict=table_dict)

        # plot 
        df_umi = pd.read_csv(self.count_file, sep='\t', index_col=None)
        filter_cells = pd.read_csv(f'{self.outs}/filtered_contig_annotations.csv', sep=',', index_col=None)
        barcode_dic_df = pd.read_csv(self.barcode_dic, sep='\t', index_col=None)
        
        filter_cells['10X'] = filter_cells['barcode'].apply(lambda x: x.strip('-1'))
        
        df = pd.merge(filter_cells, barcode_dic_df, on='10X', how='left')
        df.to_csv(f'{self.cwd_path}/{self.outdir}/contigs.csv', sep=',', index=False)
        sgr_cbs = set(df['sgr'].tolist())
        cell_nums = len(sgr_cbs)
        
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in sgr_cbs else 'UB')
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi.to_csv(self.count_file, sep='\t', index=False)
        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(self.count_file))
        
        # stat file
        assemble_summary = []

        fl_pro_pair_df = pd.DataFrame(filter_cells[(filter_cells['full_length']==True)&(filter_cells['productive']==True)].drop_duplicates(['barcode', 'chain']).barcode.value_counts())
        fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']==2]
        sum_dict = pd.read_csv(f'{self.outs}/metrics_summary.csv', sep=',', index_col=None)
        sum_dict = sum_dict.T.to_dict()
        
        read_count = int(sum_dict[0]["Number of Read Pairs"].replace(',', ''))
        
        assemble_summary.append({
            'item': 'Estimated Number of Cells',
            'count': cell_nums,
            'total_count': np.nan
        })
        assemble_summary.append({
            'item': 'Reads Mapped to Any V(D)J Gene',
            'count': int(read_count * (float(sum_dict[0]['Reads Mapped to Any V(D)J Gene'].strip('%'))/100)), 
            'total_count': read_count
        })
        assemble_summary.append({
            'item': 'Reads Mapped to TRA',
            'count': int(read_count * (float(sum_dict[0]['Reads Mapped to TRA'].strip('%'))/100)),
            'total_count': read_count,
        })
        assemble_summary.append({
            'item': 'Reads Mapped to TRB',
            'count': int(read_count * (float(sum_dict[0]['Reads Mapped to TRB'].strip('%'))/100)),
            'total_count': read_count
        })
        assemble_summary.append({
            'item': 'Mean Read Pairs per Cell',
            'count': int(sum_dict[0]['Mean Read Pairs per Cell'].replace(',', '')),
            'total_count': np.nan
        })
        assemble_summary.append({
            'item': 'Mean Used Read Pairs per Cell',
            'count': int(str(sum_dict[0]['Mean Used Read Pairs per Cell']).replace(',', '')),
            'total_count': np.nan
        })
        assemble_summary.append({
            'item': 'Median TRA UMIs per Cell',
            'count': int(float(sum_dict[0]['Median TRA UMIs per Cell'])), 
            'total_count': np.nan
        })
        assemble_summary.append({
            'item': 'Median TRB UMIs per Cell',
            'count': int(float(sum_dict[0]['Median TRB UMIs per Cell'])), 
            'total_count': np.nan
        })
        assemble_summary.append({
            'item': 'Cells With Productive V-J Spanning Pair',
            'count': fl_pro_pair_df.shape[0],
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With TRA Contig',
            'count': len(set(filter_cells[filter_cells['chain']=='TRA'].barcode.tolist())),
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With TRB Contig',
            'count': len(set(filter_cells[filter_cells['chain']=='TRB'].barcode.tolist())),
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With CDR3-annotated TRA Contig',
            'count': len(set(filter_cells[(filter_cells['cdr3']!='None')&(filter_cells['chain']=='TRA')].barcode.tolist())),
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With CDR3-annotated TRB Contig',
            'count': len(set(filter_cells[(filter_cells['cdr3']!='None')&(filter_cells['chain']=='TRB')].barcode.tolist())),
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With V-J Spanning TRA Contig',
            'count': len(set(filter_cells[(filter_cells['full_length']==True)&(filter_cells['chain']=='TRA')].barcode.tolist())),
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With V-J Spanning TRB Contig',
            'count': len(set(filter_cells[(filter_cells['full_length']==True)&(filter_cells['chain']=='TRB')].barcode.tolist())),
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With Productive TRA Contig',
            'count': len(set(filter_cells[(filter_cells['productive']==True)&(filter_cells['chain']=='TRA')].barcode.tolist())),
            'total_count': cell_nums
        })
        assemble_summary.append({
            'item': 'Cells With Productive TRB Contig',
            'count': len(set(filter_cells[(filter_cells['productive']==True)&(filter_cells['chain']=='TRB')].barcode.tolist())),
            'total_count': cell_nums
        })
        
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(assemble_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)
        
    def run(self):
        
        if not os.path.exists(self.outs):
            self.run_assemble()
        self.gen_stat_file()
        self.clean_up()
        
        
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    parser.add_argument('--species', help='species', choices=['hs','mmu'], required=True)
    parser.add_argument('--soft', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'], 
        default='4.0.0')
    parser.add_argument('--mem', help='memory (G)', default=10)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir', required=True)
        parser.add_argument('--count_file', help='count umi types for each barcode', required=True)
        parser.add_argument('--barcode_dic', help='10X barcode correspond sgr barcode', required=True)
    return parser
