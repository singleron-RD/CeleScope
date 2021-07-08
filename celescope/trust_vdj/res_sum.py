import pandas as pd
from celescope.tools.step import Step, s_common
from celescope.tools import utils
import numpy as np
import pysam
import os
from collections import defaultdict
from Bio.Seq import Seq

DIR = '/SGRNJ03/randd/zhouxin/software/TRUST4/'


class Res_sum(Step):
    """
    Features

    - Calculate clonetypes.

    Output
    - `06.res_sum/clonetypes.tsv` Record each clonetype and its frequent.
    - `06.res_sum/all_{type}.csv` Containing detailed information for each barcode.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample
        self.Seqtype = args.Seqtype
        self.filter_rep = args.filter_rep
        self.fa = args.fa

        if self.Seqtype == 'TCR':
            self.string = 't'
            self.chain = ['TRA', 'TRB']
            self.paired_groups = ['TRA_TRB']
        elif self.Seqtype == 'BCR':
            self.string = 'b'
            self.chain = ['IGH', 'IGL', 'IGK']
            self.paired_groups = ['IGH_IGL', 'IGH_IGK']
            
        self.match_bool = True
        if (not args.match_dir) or (args.match_dir == "None"):
            self.match_bool = False
        if self.match_bool:
            self.match_cell_barcodes, _match_cell_number = utils.read_barcode_file(
                args.match_dir)
 

    @utils.add_log
    def gen_stat_file(self):
        res_sum_summary = []
        data = pd.read_csv(self.filter_rep, sep=',')
        pro_data = data[data['productive'] == True]
        pro_contigs = set(pro_data['contig_id'].tolist())
        barcodes = set(data['barcode'].tolist())
        total_count = len(barcodes)
        res = pd.DataFrame()

        for c in self.chain:
            tmp = pro_data[pro_data['chain']==c]

            tmp = tmp.set_index('barcode')
            tmp = tmp.rename(columns={'chain': f'{c}_chain'})

            res = pd.concat([res, tmp], axis=1, join='outer', sort=False).fillna('None')
            
        for pg in self.paired_groups:
            attrs = pg.split('_')
            chain1 = attrs[0]
            chain2 = attrs[1]
            tmp = res[(res[f'{chain1}_chain']!='None') & (res[f'{chain2}_chain']!='None')]
            item = f'Number of cells with productive paired {chain1} and {chain2}'
            count = int(tmp.shape[0])
            res_sum_summary.append({
                'item': item,
                'count': count,
                'total_count': total_count
            })
        
        
        with pysam.FastxFile(self.fa, 'r') as fa:
            dic = defaultdict(set)

            for entry in fa:
                comment = entry.comment
                name = entry.name
                attrs = comment.split(' ')
                V_gene = attrs[2]
                D_gene = attrs[3]
                J_gene = attrs[4]
                cdr3 = attrs[7]

                if V_gene != '*':
                    chain = V_gene[:3]
                elif D_gene != '*':
                    chain = D_gene[:3]
                elif J_gene != '*':
                    chain = J_gene[:3]
                else:
                    chain = 'None'
                if chain != 'None':
                    for c in self.chain:
                        if c == chain:
                            dic[f'Cells With {c} Contig'].add(name.split('_')[0])
                            if cdr3.split('=')[1] != 'null':
                                dic[f'Cells With CDR3 {c} Contig'].add(name.split('_')[0])
                            if (V_gene != '*') and (J_gene != '*'):
                                dic[f'Cells With V-J spanning {c} Contig'].add(name.split('_')[0])
                            if name in pro_contigs:
                                dic[f'Cells With productive {c} Contig'].add(name.split('_')[0])
                        else:
                            continue
                else:
                    continue
                        
            for item in list(dic.keys()):
                res_sum_summary.append({
                    'item': item,
                    'count': len(dic[item]), 
                    'total_count': total_count
                })
                                    
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(res_sum_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)          

        if self.match_bool:
            vdj_barcodes = data['barcode'].tolist()
            rna_barcodes = []
            for barcode in self.match_cell_barcodes:
                barcode = Seq(barcode)
                barcode = barcode.reverse_complement()
                barcode = str(barcode)
                rna_barcodes.append(barcode)
                
            matched_barcodes = list(set(vdj_barcodes).intersection(set(rna_barcodes)))
            res_sum_summary.append({
                'item': 'Number of matched barcodes with scRNA-seq results',
                'count': len(matched_barcodes),
                'total_count': np.nan
            })
            
            df_bc = pd.DataFrame(matched_barcodes, columns=['barcode'])
            df_match = pd.merge(df_bc, data, on='barcode', how='inner')
            if not os.path.exists(f'{self.outdir}/matched_res/'):
                os.makedirs(f'{self.outdir}/matched_res/')
            df_match.to_csv(f'{self.outdir}/matched_res/filter_contig.csv', sep=',', index=False)
            
            
    @utils.add_log
    def run(self):

        self.gen_stat_file()

        self.clean_up()


@utils.add_log
def res_sum(args):
    step_name = 'res_sum'
    res_sum_obj = Res_sum(args, step_name)
    res_sum_obj.run()


def get_opts_res_sum(parser, sub_program):
    parser.add_argument('--Seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--filter_rep', help='filtered assemble report without imputation', required=True)
        parser.add_argument('--fa', help='assembled fasta file', required=True)
        parser.add_argument('--match_dir', help='scRNA-seq results directory', default=None)







    
