import numpy as np
import pandas as pd
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common


class Match(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        
        self.match_dir = args.match_dir
        self.contig_df = args.contig_df
        
        # match_dir
        self.match_bool = True
        if (not args.match_dir) or (args.match_dir == "None"):
            self.match_bool = False
        if self.match_bool:
            self.match_cell_barcodes, _match_cell_number = utils.read_barcode_file(
                args.match_dir)
    
    @utils.add_log     
    def run_match(self):
        df_sgr = pd.DataFrame(self.match_cell_barcodes, columns=['sgr'])
        contig_df = pd.read_csv(self.contig_df, sep=',', index_col=None)
        contig_df['sgr'] = contig_df['sgr'].apply(lambda x: str(Seq(x).reverse_complement ()))
        df = pd.merge(contig_df, df_sgr, on='sgr', how='inner')
        df.to_csv(f'{self.outdir}/match_contigs.csv', sep=',')
        return df
    
    @utils.add_log
    def gen_stat_file(self):
        df = self.run_match()
        cell_nums = len(set(df['sgr'].tolist()))
        fl_pro_pair_df = pd.DataFrame(df[(df['full_length']==True)&(df['productive']==True)].drop_duplicates(['barcode', 'chain']).barcode.value_counts())
        fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']==2]
        
        match_summary = []
        match_summary.append({
            'item': 'Cells consistent with scRNA-seq analysis',
            'count': len(set(df.sgr.tolist())),
            'total_count': np.nan
        })
        match_summary.append({
            'item': 'Median TRA UMIs per cell',
            'count': int(np.median(df[df['chain']=='TRA']['umis'].tolist())),
            'total_count': np.nan
        })
        match_summary.append({
            'item': 'Median TRB UMIs per cell',
            'count': int(np.median(df[df['chain']=='TRB']['umis'].tolist())),
            'total_count': np.nan
        })
        match_summary.append({
            'item': 'Cells With Productive V-J Spanning Pair',
            'count': fl_pro_pair_df.shape[0],
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With TRA Contig',
            'count': len(set(df[df['chain']=='TRA'].barcode.tolist())),
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With TRB Contig',
            'count': len(set(df[df['chain']=='TRB'].barcode.tolist())),
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With CDR3-annotated TRA Contig',
            'count': len(set(df[(df['cdr3']!='None')&(df['chain']=='TRA')].barcode.tolist())),
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With CDR3-annotated TRB Contig',
            'count': len(set(df[(df['cdr3']!='None')&(df['chain']=='TRB')].barcode.tolist())),
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With V-J Spanning TRA Contig',
            'count': len(set(df[(df['full_length']==True)&(df['chain']=='TRA')].barcode.tolist())),
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With V-J Spanning TRB Contig',
            'count': len(set(df[(df['full_length']==True)&(df['chain']=='TRB')].barcode.tolist())),
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With Productive TRA Contig',
            'count': len(set(df[(df['productive']==True)&(df['chain']=='TRA')].barcode.tolist())),
            'total_count': cell_nums
        })
        match_summary.append({
            'item': 'Cells With Productive TRB Contig',
            'count': len(set(df[(df['productive']==True)&(df['chain']=='TRB')].barcode.tolist())),
            'total_count': cell_nums
        })
        
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(match_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)
        
    def run(self):
        self.gen_stat_file()
        self.clean_up()
        
def match(args):
    step_name = 'match'
    match_obj = Match(args, step_name)
    match_obj.run()
    
def get_opts_match(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq analysis directory', required=True)
        parser.add_argument('--contig_df', help='contig dataframes with sgr barcodes', required=True)

        

        
        