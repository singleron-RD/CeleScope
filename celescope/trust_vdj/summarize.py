from collections import defaultdict

import pandas as pd
import numpy as np
import pysam
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.trust_vdj.__init__ import CHAIN, PAIRED_CHAIN
from celescope.tools.cellranger3 import get_plot_elements


def get_vj_annot(df, chains, pairs):
    fl_pro_pair_df = pd.DataFrame(df[df['productive']==True].barcode.value_counts())
    fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode']>=2]
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
            'count': len(set(df[(df['chain']==c)&(df['productive']==True)].barcode.tolist())),
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
        self.reads_assignment = args.reads_assignment
        self.fq2 = args.fq2
        self.assembled_fa = args.assembled_fa

        self.chains = CHAIN[self.seqtype]
        self.paired_groups = PAIRED_CHAIN[self.seqtype]
        

        # common variables

    
    @utils.add_log
    def process(self):
        df = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/{self.sample}_contig.csv', sep='\t', header=None)
        df.columns = ['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id', 'raw_consensus_id']

        df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
        df['productive'] = df['cdr3'].apply(lambda x: False if x=='None' else True)

        df.to_csv(f'{self.outdir}/{self.sample}_contig.csv', sep=',', index=False)
        # df = df.sort_values(by='umis', ascending=False)
        if self.seqtype == 'BCR':
            igh = df[df['chain']=='IGH']
            temp = df[(df['chain']=='IGK') | (df['chain']=='IGL')]
            temp = temp.sort_values(by='umis', ascending=False)
            temp = temp.drop_duplicates(['barcode'])
            df_for_clono = pd.concat([igh, temp], ignore_index=True)
            df_for_clono = df_for_clono.sort_values(by='umis', ascending=False)
            df_for_clono = df_for_clono.drop_duplicates(['barcode', 'chain'])
        else:
            df = df[(df['chain']=='TRA') | (df['chain']=='TRB')]
            df = df.sort_values(by='umis', ascending=False)
            df_for_clono = df.drop_duplicates(['barcode', 'chain'])
        # gen clonotypes table
        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes = set(df_for_clono_pro['barcode'].tolist())
        total_cells =len(cell_barcodes)

        summarize_summary = []
        summarize_summary.append({
            'item': 'Estimated Number of Cells',
            'count': total_cells,
            'total_count': np.nan
        })

        df_for_clono_pro['chain_cdr3aa'] = df_for_clono_pro[['chain', 'cdr3']].apply(':'.join, axis=1)   

        cbs = set(df_for_clono_pro['barcode'].tolist())
        clonotypes = open(f'{self.outdir}/clonotypes.csv', 'w')
        clonotypes.write('barcode\tcdr3s_aa\n')
        for cb in cbs:
            temp = df_for_clono_pro[df_for_clono_pro['barcode']==cb]
            temp = temp.sort_values(by='chain', ascending=True)
            chain_list = temp['chain_cdr3aa'].tolist()
            chain_str = ';'.join(chain_list)
            clonotypes.write(f'{cb}\t{chain_str}\n')
        clonotypes.close() 

        df_clonotypes = pd.read_csv(f'{self.outdir}/clonotypes.csv', sep='\t', index_col=None)
        df_clonotypes = df_clonotypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
        df_clonotypes = df_clonotypes.rename(columns={'barcode': 'frequency'})
        sum_f = df_clonotypes['frequency'].sum()
        df_clonotypes['proportion'] = df_clonotypes['frequency'].apply(lambda x: x/sum_f)
        df_clonotypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonotypes.shape[0]+1)]
        df_clonotypes = df_clonotypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
        df_clonotypes = df_clonotypes.sort_values(by='frequency', ascending=False)
        df_clonotypes.to_csv(f'{self.outdir}/clonotypes.csv', sep=',', index=False) 

        df_clonotypes['ClonotypeID'] = df_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_clonotypes['Frequency'] = df_clonotypes['frequency']
        df_clonotypes['Proportion'] = df_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        df_clonotypes['CDR3_aa'] = df_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', df_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']])
        self.add_data_item(table_dict=table_dict)

        # reads summary
        read_count = 0
        umi_dict = defaultdict(set)
        umi_count = defaultdict()
        with pysam.FastxFile(self.fq2) as fq:
            for read in fq:
                read_count+=1
                cb = read.name.split('_')[0]
                umi = read.name.split('_')[1]
                umi_dict[cb].add(umi)
        for cb in umi_dict:
            umi_count[cb] = len(umi_dict[cb])
        df_umi = pd.DataFrame.from_dict(umi_count, orient='index', columns=['UMI'])
        df_umi['barcode'] = df_umi.index
        df_umi = df_umi.reset_index(drop=True)
        df_umi = df_umi.reindex(columns=['barcode', 'UMI'])
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in cell_barcodes else 'UB')
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)

        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))
        

        summarize_summary.append({
            'item': 'Mean Read Pairs per Cell',
            'count': int(read_count/total_cells),
            'total_count': np.nan
        })
        with pysam.FastaFile(self.assembled_fa) as fa:
            summarize_summary.append({
                'item': 'Mean Used Read Pairs per Cell',
                'count': int(fa.nreferences/total_cells), 
                'total_count': np.nan
            })
            summarize_summary.append({
                'item': 'Fraction of Reads in Cells',
                'count': fa.nreferences,
                'total_count': read_count
            })
        
        for c in self.chains:
            temp_df = df_for_clono_pro[df_for_clono_pro['chain']==c]
            summarize_summary.append({
                'item': f'Median {c} UMIs per Cell',
                'count': int(temp_df['umis'].median()),
                'total_count': np.nan
            })
        
        annotation_summary = get_vj_annot(df_for_clono, self.chains, self.paired_groups)
        summarize_summary = summarize_summary + annotation_summary
        # gen stat file          
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(summarize_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)     

        
    @utils.add_log
    def run(self):
        self.process()
        self.clean_up()


@utils.add_log
def summarize(args):
    step_name = f'{args.seqtype}_summarize'
    summarize_obj = Summarize(args, step_name)
    summarize_obj.run()


def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--reads_assignment', help='File records reads assigned to contigs.', required=True)
        parser.add_argument('--fq2', help='Cutadapt R2 reads.', required=True)
        parser.add_argument('--assembled_fa', help='Read used for assembly', required=True)








    
