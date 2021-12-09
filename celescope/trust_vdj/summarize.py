from collections import defaultdict

import pandas as pd
import numpy as np
import pysam
import copy
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.trust_vdj.__init__ import CHAIN, PAIRED_CHAIN
from celescope.tools.cellranger3 import get_plot_elements
import sys

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

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

    - TCR/BCR full length assembly results.

    Output
    - `04.summarize/clonetypes.tsv` High-level descriptions of each clonotype.
    - `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
    - `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
    - `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig_annotations.csv.
    - `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
    - `04.summarize/{sample}_chain_filtered_contig.csv`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.csv.
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
        
        # use for cell filtering  default = 4
        self.min_read_count = args.min_read_count

        # common variables

    
    @utils.add_log
    def process(self):
        df = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/{self.sample}_contig.csv', sep='\t', header=None)
        df.columns = ['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis']

        df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
        df['productive'] = df['cdr3'].apply(lambda x: False if x=='None' else True)

        # all contig.csv
        df_revcompl = copy.deepcopy(df)
        df_revcompl['barcode'] = df_revcompl['barcode'].apply(lambda x: reversed_compl(x))
        df_revcompl['contig_id'] = df_revcompl['contig_id'].apply(lambda x: reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])
        
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
            df_for_clono = df[(df['chain']=='TRA') | (df['chain']=='TRB')]
            df_for_clono = df_for_clono.sort_values(by='umis', ascending=False)
            df_for_clono = df_for_clono.drop_duplicates(['barcode', 'chain'])

        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes = set(df_for_clono_pro['barcode'].tolist())
        total_cells_all =len(cell_barcodes) # cell number before filter
        
        '''
        filter cell by cell type from trust barcode filtered report
        filter cell by read count from trust report 
        '''
        filterbc = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/barcoderepfl.tsv',sep='\t')
        trust_rep = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/report.out', sep='\t')
        filterbc = filterbc.rename(columns = {'#barcode':'barcode'})
        filterbc = filterbc[(filterbc['chain1']!='*') | (filterbc['chain2']!='*')]
        filterbc = set(filterbc['barcode'].tolist())

        trust_rep = trust_rep[trust_rep['cid_full_length'] >= 1]
        trust_rep = trust_rep[trust_rep['CDR3aa'] != 'out_of_frame']
        trust_rep = trust_rep.rename(columns = {'cid':'barcode', '#count':'count'})
        trust_rep['barcode'] = trust_rep['barcode'].apply(lambda x:x.split('_')[0])
        trust_rep = trust_rep.sort_values(by = 'count', ascending = False)
        if self.min_read_count == "auto":
            trust_rep = trust_rep[trust_rep['count'] >= 4]
        else:
            min_read_count = int(self.min_read_count)
            trust_rep = trust_rep[trust_rep['count'] >= min_read_count]
        trust_rep = set(trust_rep['barcode'].tolist())

        df_for_clono = df_for_clono[df_for_clono['barcode'].isin(filterbc)]
        df_for_clono = df_for_clono[df_for_clono['barcode'].isin(trust_rep)]


        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes = set(df_for_clono_pro['barcode'].tolist())
        total_cells =len(cell_barcodes)
        
        # filter contig.csv
        df_filter_contig = copy.deepcopy(df)
        df_filter_contig = df_filter_contig[df_filter_contig['barcode'].isin(cell_barcodes)]
        df_filter_contig['barcode'] = df_filter_contig['barcode'].apply(lambda x: reversed_compl(x))
        df_filter_contig['contig_id'] = df_filter_contig['contig_id'].apply(lambda x: reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])

        # filter contig.fasta
        all_contig_fasta = f'{self.outdir}/{self.sample}_all_contig.fasta'
        filter_contig_fasta = f'{self.outdir}/{self.sample}_filtered_contig.fasta'
        filter_contig_fasta = open(filter_contig_fasta,'w')
        with pysam.FastxFile(all_contig_fasta) as fa:
            for read in fa:
                name = read.name
                barcode = name.split('_')[0]
                sequence = read.sequence
                if reversed_compl(barcode) in cell_barcodes:
                    filter_contig_fasta.write('>' + name + '\n' + sequence + '\n')
        filter_contig_fasta.close()


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
        contig_with_clonotype = copy.deepcopy(df_clonotypes)

        df_clonotypes = df_clonotypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
        df_clonotypes = df_clonotypes.rename(columns={'barcode': 'frequency'})
        sum_f = df_clonotypes['frequency'].sum()
        df_clonotypes['proportion'] = df_clonotypes['frequency'].apply(lambda x: x/sum_f)
        df_clonotypes = df_clonotypes.sort_values(by='frequency', ascending=False)
        df_clonotypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonotypes.shape[0]+1)]
        df_clonotypes = df_clonotypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
        df_clonotypes.to_csv(f'{self.outdir}/clonotypes.csv', sep=',', index=False) 
        used_for_merge = df_clonotypes[['cdr3s_aa','clonotype_id']]

        df_clonotypes['ClonotypeID'] = df_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_clonotypes['Frequency'] = df_clonotypes['frequency']
        df_clonotypes['Proportion'] = df_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        df_clonotypes['CDR3_aa'] = df_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', df_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']])
        self.add_data_item(table_dict=table_dict)

        # add clonotype_id in contig.csv file
        df_merge = pd.merge(used_for_merge, contig_with_clonotype, on='cdr3s_aa', how='outer')
        df_merge = df_merge[['barcode','clonotype_id']]
        df_merge['barcode'] = df_merge['barcode'].apply(lambda x: reversed_compl(x))
        df_revcompl = pd.merge(df_merge, df_revcompl, on='barcode',how='outer')
        df_filter_contig = pd.merge(df_merge, df_filter_contig, on='barcode',how='outer')
        df_revcompl.fillna('None',inplace = True)
        df_filter_contig.fillna('None',inplace = True)
        df_revcompl = df_revcompl[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_filter_contig = df_filter_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_revcompl.to_csv(f'{self.outdir}/{self.sample}_all_contig.csv', sep=',', index=False)
        df_filter_contig.to_csv(f'{self.outdir}/{self.sample}_filtered_contig.csv', sep=',', index=False)

        # filter chains based on umis
        df_filter_contig.sort_values(by=['barcode','umis'],ascending=[True,False],inplace=True)
        df_filter_contig.reset_index(drop=True,inplace=True)
        chain_filter_dict = df_filter_contig.groupby('barcode')['contig_id'].apply(lambda x:x.tolist()).to_dict()
        chain_filter_set = set()
        for key,val in chain_filter_dict.items():
            if len(val)>2:
                for contig_id in val[:2]:
                    chain_filter_set.add(contig_id)
            else:
                chain_filter_set.add(val[0])
        df_filter_contig = df_filter_contig[df_filter_contig['contig_id'].isin(chain_filter_set)]
        df_filter_contig.to_csv(f'{self.outdir}/{self.sample}_chain_filtered_contig.csv', sep=',', index=False)

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
        df_umi['barcode'] = df_umi['barcode'].apply(lambda x: reversed_compl(x))
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)

        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))
        

        summarize_summary.append({
            'item': 'Mean Read Pairs per Cell',
            'count': int(read_count/total_cells_all),
            'total_count': np.nan
        })
        with pysam.FastaFile(self.assembled_fa) as fa:
            summarize_summary.append({
                'item': 'Mean Used Read Pairs per Cell',
                'count': int(fa.nreferences/total_cells_all), 
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
    parser.add_argument('--min_read_count', help ='filter cell by read count number, int type required', default='auto')
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--reads_assignment', help='File records reads assigned to contigs.', required=True)
        parser.add_argument('--fq2', help='Cutadapt R2 reads.', required=True)
        parser.add_argument('--assembled_fa', help='Read used for assembly', required=True)








    
