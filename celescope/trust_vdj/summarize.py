import os
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.trust_vdj.__init__ import TRUST, CHAIN, PAIRED_CHAIN


@utils.add_log
def fa_to_csv(annot_fa, assign_file, outdir, chains, average_coverage=5):
    # reads assignment 
    assignment = pd.read_csv(assign_file, sep='\t', header=None)
    contigs = open(f'{outdir}/all_contigs.csv', 'w')
    contigs.write('barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id,coverage\n')
    with pysam.FastxFile(annot_fa) as fa:
        for read in fa:
            name = read.name
            comment = read.comment
            attrs = comment.split(' ')
            barcode = name.split('_')[0]
            is_cell = 'True'
            high_confidence = 'True'
            length = attrs[0]
            if float(attrs[1]) < int(average_coverage):
                continue
            for c in chains:
                if c in comment:
                    chain = c
                    break
                else:
                    continue
            full_length = 'True'
            if not attrs[2]=='*':
                v_gene = attrs[2].split('(')[0]
            else:
                v_gene = 'None'
                full_length = 'False'
            if not attrs[3]=='*':
                d_gene = attrs[3].split('(')[0]
            else:
                d_gene = 'None'
            if not attrs[4]=='*':
                j_gene = attrs[4].split('(')[0]
            else:
                j_gene = 'None'
                full_length = 'False'
            if not attrs[5]=='*':
                c_gene = attrs[5].split('(')[0]
            else:
                c_gene = 'None'

            cdr3 = attrs[8].split('=')[1]
            cdr3_aa = 'None'
            productive = 'False'
            if not cdr3 == 'null':
                cdr3_aa = str(Seq(cdr3).translate())
                if (int(len(cdr3)) % 3 == 0) and (not '*' in cdr3_aa):
                    productive = 'True'
            temp = assignment[assignment[1]==name]
            read_list = [i for i in temp[0].tolist() if i.split('_')[0] in name]
            reads = str(len(read_list))
            umis = str(len(set([i.split("_")[1] for i in read_list])))
            raw_consensus_id = 'None'
            raw_clonotype_id = 'None'
                
            string = ','.join([barcode, is_cell, name, high_confidence, length, chain, v_gene, d_gene, j_gene, c_gene, full_length, productive, cdr3_aa, cdr3, reads, umis, raw_clonotype_id, raw_consensus_id, attrs[1]])
            contigs.write(f'{string}\n')

    contigs.close()
    df_all = pd.read_csv(f'{outdir}/all_contig.csv', sep=',')
    df_all = df_all.sort_values(by='coverage', ascending=False)
    df_all = df_all.drop_duplicates(['barcode', 'chain'])

    return df_all


def get_vj_annot(df, chains, pairs):
    fl_pro_pair_df = pd.DataFrame(df[(df['full_length']==True)&(df['productive']==True)].drop_duplicates(['barcode', 'chain']).barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']==2]
    l = []
    cell_nums = len(set(df['barcode'].tolist()))
    l.append({
        'item': 'Cells With Productive V-J Spanning Pair',
        'count': fl_pro_pair_df.shape[0],
        'total_count': cell_nums
    })
    for p in pairs:
        chain1 = p.split('_')[0]
        chain2 = p.split('_')[1]
        cbs1 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain1)].barcode.tolist())
        cbs2 = set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==chain2)].barcode.tolist())
        paired_cbs = len(cbs1.intersection(cbs2))
        l.append({
            'item': f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
            'count': paired_cbs,
            'total_count': cell_nums
        })
    for c in chains:
        l.append({
            'item': f'Cells With {c} Contig',
            'count': len(set(df[df['chain']==c].barcode.tolist())),
            'total_count': cell_nums
        })
        l.append({
            'item': f'Cells With CDR3-annotated {c} Contig',
            'count': len(set(df[(df['chain']==c)&(df['cdr3']!='None')].barcode.tolist())),
            'total_count': cell_nums
        })
        l.append({
            'item': f'Cells With V-J Spanning {c} Contig',
            'count': len(set(df[(df['full_length']==True)&(df['chain']==c)].barcode.tolist())),
            'total_count': cell_nums
        })
        l.append({
            'item': f'Cells With Productive {c} Contig',
            'count': len(set(df[(df['full_length']==True)&(df['productive']==True)&(df['chain']==c)].barcode.tolist())),
            'total_count': cell_nums
        })

    return l

class Summarize(Step):
    """
    Features

    - Calculate clonetypes.

    Output
    - `06.summarize/clonetypes.tsv` Record each clonetype and its frequent.
    - `06.summarize/all_{type}.csv` Containing detailed information for each barcode.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample
        self.seqtype = args.seqtype
        self.all_rep = args.all_rep
        self.fa = args.fa
        self.annot_fa = args.annot_fa

        self.chains = CHAIN[self.seqtype]
        self.paired_groups = PAIRED_CHAIN[self.seqtype]

        # common variables
            
    @utils.add_log
    def get_full_len_assembly(self):
        cmd = (
            f'perl {TRUST}/GetFullLengthAssembly.pl {self.annot_fa} > {self.outdir}/all_contigs.fasta'
        )
        Summarize.get_full_len_assembly.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    
    @utils.add_log
    def filter_table(self):
        df = pd.read_csv(f'{self.outdir}/all_contigs.csv', sep=',')
        df_sort = df.sort_values(by='reads', ascending=False)
        df_filter = df_sort.drop_duplicates(['barcode', 'chain'])
        df_filter.to_csv(f'{self.outdir}/filtered_contigs.csv')


    @utils.add_log
    def gen_stat_file(self):
        summarize_summary = []
        data = pd.read_csv(self.all_rep, sep=',')
        pro_data = data[(data['full_length']==True) & (data['productive']==True)]
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
            item = f'Number of cells with V-J spanning productive paired {chain1} and {chain2}'
            count = int(tmp.shape[0])
            summarize_summary.append({
                'item': item,
                'count': count,
                'total_count': total_count
            })
        
        for c in self.chain:
            dic = {}
            tmp = data[data['chain'] == c]
            dic[f'Cells With {c} Contig'] = tmp.shape[0]
            dic[f'Cells With V-J spanning {c} Contig'] = tmp[tmp['full_length']==True].shape[0]
            dic[f'Cells With productive {c} Contig'] = tmp[(tmp['productive']==True) & (tmp['full_length']==True)].shape[0]
                        
            for item in list(dic.keys()):
                summarize_summary.append({
                    'item': item,
                    'count': int(dic[item]), 
                    'total_count': total_count
                })
                                    
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(summarize_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)          

            
    @utils.add_log
    def run(self):

        self.gen_stat_file()

        self.clean_up()


@utils.add_log
def summarize(args):
    step_name = 'summarize'
    summarize_obj = Summarize(args, step_name)
    summarize_obj.run()


def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)

    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--all_rep', help='filtered assemble report without imputation', required=True)
        parser.add_argument('--annot_fa', help='assembled fasta contig file.', required=True)








    
