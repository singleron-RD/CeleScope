from collections import defaultdict
from collections import Counter
import pandas as pd
import numpy as np
import pysam
import copy
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.trust_vdj.__init__ import CHAIN, PAIRED_CHAIN
from celescope.tools.cellranger3 import get_plot_elements
from celescope.trust_vdj import trust_utils as tr

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
    - `04.summarize/{sample}_chain_filtered_contig.fasta`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.fasta.
    - `04.summarize/{sample}_one_chain_contig.csv`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.csv.
    - `04.summarize/{sample}_one_chain_contig.fasta`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.fasta.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample
        self.seqtype = args.seqtype
        self.reads_assignment = args.reads_assignment
        self.fq2 = args.fq2
        self.assembled_fa = args.assembled_fa

        self.chains, self.paired_groups = self.parse_seqtype(self.seqtype)
        self.min_read_count = args.min_read_count
        self.summarize_summary = []

    @staticmethod
    def reversed_compl(seq):
        return str(Seq(seq).reverse_complement())
    
    @staticmethod
    def parse_seqtype(seqtype):
        return CHAIN[seqtype], PAIRED_CHAIN[seqtype]
    
    @staticmethod
    def filter_cell(df, seqtype, filterbc_rep, trust_rep, min_read_count):
        if seqtype == 'BCR':
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

        filterbc_rep = filterbc_rep.rename(columns = {'#barcode':'barcode'})
        filterbc = filterbc_rep[(filterbc_rep['chain1']!='*') | (filterbc_rep['chain2']!='*')]
        filterbc = set(filterbc['barcode'].tolist())

        trust_rep = trust_rep[trust_rep['cid_full_length'] >= 1]
        trust_rep = trust_rep[trust_rep['CDR3aa'] != 'out_of_frame']
        trust_rep = trust_rep.rename(columns = {'cid':'barcode', '#count':'count'})
        trust_rep['barcode'] = trust_rep['barcode'].apply(lambda x:x.split('_')[0])
        trust_rep = trust_rep.sort_values(by = 'count', ascending = False)
        if min_read_count == "auto":
            trust_rep = trust_rep[trust_rep['count'] >= 4]
        else:
            min_read_count = int(min_read_count)
            trust_rep = trust_rep[trust_rep['count'] >= min_read_count]
        trust_rep_bc = set(trust_rep['barcode'].tolist())

        df_for_clono = df_for_clono[df_for_clono['barcode'].isin(filterbc)]
        df_for_clono = df_for_clono[df_for_clono['barcode'].isin(trust_rep_bc)]

        return df_for_clono
    
    @staticmethod
    def init_contig(df, cell_barcodes):
        # all contig.csv
        df_all_contig = copy.deepcopy(df)
        df_all_contig['barcode'] = df_all_contig['barcode'].apply(lambda x: Summarize.reversed_compl(x))
        df_all_contig['contig_id'] = df_all_contig['contig_id'].apply(lambda x: Summarize.reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])
        # filter contig.csv
        df_filter_contig = copy.deepcopy(df)
        df_filter_contig = df_filter_contig[df_filter_contig['barcode'].isin(cell_barcodes)]
        df_filter_contig['barcode'] = df_filter_contig['barcode'].apply(lambda x: Summarize.reversed_compl(x))
        df_filter_contig['contig_id'] = df_filter_contig['contig_id'].apply(lambda x: Summarize.reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])

        return df_all_contig, df_filter_contig
    
    @staticmethod
    def filter_fasta(cell_barcodes, all_contig_fasta, filter_contig_fasta):
        with pysam.FastxFile(all_contig_fasta) as fa:
            for read in fa:
                name = read.name
                barcode = name.split('_')[0]
                sequence = read.sequence
                if Summarize.reversed_compl(barcode) in cell_barcodes:
                    filter_contig_fasta.write('>' + name + '\n' + sequence + '\n')
        filter_contig_fasta.close()

    @staticmethod
    def parse_clonotypes(df_for_clono_pro, outdir):
        df_for_clono_pro['chain_cdr3aa'] = df_for_clono_pro[['chain', 'cdr3']].apply(':'.join, axis=1)
        df_for_clono_pro['chain_cdr3nt'] = df_for_clono_pro[['chain', 'cdr3_nt']].apply(':'.join, axis=1)

        cbs = set(df_for_clono_pro['barcode'].tolist())
        clonotypes = open(f'{outdir}/clonotypes.csv', 'w')
        clonotypes.write('barcode\tcdr3s_aa\tcdr3s_nt\n')
        for cb in cbs:
            temp = df_for_clono_pro[df_for_clono_pro['barcode']==cb]
            temp = temp.sort_values(by='chain', ascending=True)
            chain_list = temp['chain_cdr3aa'].tolist()
            chain_str = ';'.join(chain_list)
            ntchain_list = temp['chain_cdr3nt'].tolist()
            ntchain_str = ';'.join(ntchain_list)
            clonotypes.write(f'{cb}\t{chain_str}\t{ntchain_str}\n')
        clonotypes.close() 

        df_clonotypes = pd.read_csv(f'{outdir}/clonotypes.csv', sep='\t', index_col=None)
        df_dict = df_clonotypes[["cdr3s_nt", "cdr3s_aa"]].set_index("cdr3s_nt").to_dict(orient='dict')['cdr3s_aa']
        contig_with_clonotype = copy.deepcopy(df_clonotypes)

        df_clonotypes = df_clonotypes.groupby('cdr3s_nt', as_index=False).agg({'barcode': 'count'})
        df_clonotypes = df_clonotypes.rename(columns={'barcode': 'frequency'})
        sum_f = df_clonotypes['frequency'].sum()
        df_clonotypes['proportion'] = df_clonotypes['frequency'].apply(lambda x: x/sum_f)
        df_clonotypes = df_clonotypes.sort_values(by='frequency', ascending=False)
        df_clonotypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonotypes.shape[0]+1)]
        df_clonotypes['cdr3s_aa'] = df_clonotypes['cdr3s_nt'].apply(lambda x:df_dict[x])
        df_clonotypes = df_clonotypes.reindex(columns=['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt'])
        df_clonotypes.to_csv(f'{outdir}/clonotypes.csv', sep=',', index=False) 
        used_for_merge = df_clonotypes[['cdr3s_nt','clonotype_id']]

        df_clonotypes['ClonotypeID'] = df_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        df_clonotypes['Frequency'] = df_clonotypes['frequency']
        df_clonotypes['Proportion'] = df_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        df_clonotypes['CDR3_aa'] = df_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

        return df_clonotypes, contig_with_clonotype, used_for_merge

    @staticmethod
    def add_clonotypes(used_for_merge, contig_with_clonotype, df_all_contig, df_filter_contig, outdir, sample):
        # add clonotype_id in contig.csv file
        df_merge = pd.merge(used_for_merge, contig_with_clonotype, on='cdr3s_nt', how='outer')
        df_merge = df_merge[['barcode','clonotype_id']]
        df_merge['barcode'] = df_merge['barcode'].apply(lambda x: Summarize.reversed_compl(x))
        df_all_contig = pd.merge(df_merge, df_all_contig, on='barcode',how='outer')
        df_filter_contig = pd.merge(df_merge, df_filter_contig, on='barcode',how='outer')
        df_all_contig.fillna('None',inplace = True)
        df_filter_contig.fillna('None',inplace = True)
        df_all_contig = df_all_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_filter_contig = df_filter_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_all_contig.to_csv(f'{outdir}/{sample}_all_contig.csv', sep=',', index=False)
        df_filter_contig.to_csv(f'{outdir}/{sample}_filtered_contig.csv', sep=',', index=False)
        
        return df_filter_contig

    @staticmethod
    def keep_two_chains(df_filter_contig, outdir, sample):
        df_filter_contig=df_filter_contig[df_filter_contig['productive']==True]
        df_filter_contig.sort_values(by=['barcode','umis'],ascending=[True,False],inplace=True)
        df_filter_contig.reset_index(drop=True,inplace=True)
        contig_dict = df_filter_contig.groupby('barcode')['contig_id'].apply(lambda x:x.tolist()).to_dict()
        chain_dict = df_filter_contig.groupby('barcode')['chain'].apply(lambda x:x.tolist()).to_dict()
        chain_filter_set = set()
        for key,val in chain_dict.items():
            type_count=dict(Counter(val))
            for _type, _count in type_count.items():
                if _count>=2:
                    index = [_i for _i,_j in enumerate(val) if _j == _type]
                    for _index in index[:2]:
                        chain_filter_set.add(contig_dict[key][_index])
                elif _count==1:
                    index = val.index(_type)
                    chain_filter_set.add(contig_dict[key][index])
                else:
                    continue
        two_chains_contig = df_filter_contig[df_filter_contig['contig_id'].isin(chain_filter_set)]
        two_chains_contig.to_csv(f'{outdir}/{sample}_chain_filtered_contig.csv', sep=',', index=False)
        
        #fasta
        filter_contig_set = set(df_filter_contig['contig_id'].tolist())
        two_chains_fasta = f'{outdir}/{sample}_chain_filtered_contig.fasta'
        two_chains_fasta = open(two_chains_fasta,'w')
        filter_contig_fasta = f'{outdir}/{sample}_filtered_contig.fasta'
        with pysam.FastxFile(filter_contig_fasta) as fa:
            for read in fa:
                seq = read.sequence
                name = read.name
                if name in filter_contig_set:
                    two_chains_fasta.write(">"+name+"\n"+seq+"\n")
        two_chains_fasta.close()

        return two_chains_contig, two_chains_fasta
    
    @staticmethod
    def keep_one_pair(seqtype, two_chains_contig, outdir, sample):
        if seqtype == 'BCR':
            df_one_chain_heavy = two_chains_contig[ (two_chains_contig['chain'] == 'IGH') ] 
            df_one_chain_light = two_chains_contig[ (two_chains_contig['chain'] == 'IGL')|(two_chains_contig['chain'] =='IGK')]
            
            df_one_chain_heavy = df_one_chain_heavy.sort_values(by='umis', ascending=False)
            df_one_chain_light = df_one_chain_light.sort_values(by='umis', ascending=False)
            df_one_chain_heavy = df_one_chain_heavy.drop_duplicates(['barcode'])
            df_one_chain_light = df_one_chain_light.drop_duplicates(['barcode'])
            df_one_chain_merge = pd.concat([df_one_chain_heavy, df_one_chain_light], ignore_index=True)
            df_one_chain_merge = df_one_chain_merge.sort_values(by=['barcode','umis'], ascending=[True,False])
        else:
            df_one_chain = two_chains_contig[(two_chains_contig['chain']=='TRA') | (two_chains_contig['chain']=='TRB')]
            df_one_chain = df_one_chain.sort_values(by=['barcode','umis'], ascending=[True,False])
            df_one_chain_merge = df_one_chain.drop_duplicates(['barcode', 'chain'])
        df_one_chain_merge.to_csv(f'{outdir}/{sample}_one_chain_contig.csv', sep=',', index=False)
            
        # # IGH+IGK/IGL TRA+TRB pair fasta
        one_chain_contig_set = set(df_one_chain_merge['contig_id'].tolist())
        one_chain_fasta = f'{outdir}/{sample}_one_chain_contig.fasta'
        one_chain_fasta = open(one_chain_fasta,'w')
        with pysam.FastxFile(f'{outdir}/{sample}_filtered_contig.fasta') as fa:
            for read in fa:
                seq = read.sequence
                name = read.name
                if name in one_chain_contig_set:
                    one_chain_fasta.write(">"+name+"\n"+seq+"\n")
        one_chain_fasta.close()

    @utils.add_log
    def parse_contig_file(self):
        df = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/{self.sample}_contig.csv', sep='\t', header=None)
        df.columns = ['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis']
        df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
        df['productive'] = df['cdr3'].apply(lambda x: False if x=='None' else True)
    
        '''
        filter cell by cell type from trust barcode filtered report
        filter cell by read count from trust report 
        '''
        filterbc_rep = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/barcoderepfl.tsv',sep='\t')
        trust_rep = pd.read_csv(f'{self.outdir}/../03.assemble/assemble/report.out', sep='\t')
        df_for_clono = self.filter_cell(df, self.seqtype, filterbc_rep, trust_rep, self.min_read_count)

        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes = set(df_for_clono_pro['barcode'].tolist())
        total_cells =len(cell_barcodes)

        # filter contig.fasta
        all_contig_fasta = f'{self.outdir}/{self.sample}_all_contig.fasta'
        filter_contig_fasta = f'{self.outdir}/{self.sample}_filtered_contig.fasta'
        filter_contig_fasta = open(filter_contig_fasta,'w')
        self.filter_fasta(cell_barcodes, all_contig_fasta, filter_contig_fasta)

        self.summarize_summary.append({
            'item': 'Estimated Number of Cells',
            'count': total_cells,
            'total_count': np.nan
        })

        df_all_contig, df_filter_contig = self.init_contig(df, cell_barcodes)
        df_clonotypes, contig_with_clonotype, used_for_merge = self.parse_clonotypes(df_for_clono_pro,self.outdir)

        title = 'Clonetypes'
        table_dict = self.get_table(title, 'clonetypes_table', df_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']])
        self.add_data_item(table_dict=table_dict)

        # add clonotype_id in contig.csv file
        df_filter_contig = self.add_clonotypes(used_for_merge, contig_with_clonotype, df_all_contig, df_filter_contig, self.outdir, self.sample)

        # keep chains with two highest umi contigs 
        two_chains_contig, two_chains_fasta = self.keep_two_chains(df_filter_contig, self.outdir, self.sample)
        # keep one chain, IGH+IGK/IGL TRA+TRB pair
        self.keep_one_pair(self.seqtype, two_chains_contig, self.outdir, self.sample)

        return df_for_clono, df_for_clono_pro, cell_barcodes, total_cells
    
    @utils.add_log
    def gen_summary(self, df_for_clono, df_for_clono_pro, cell_barcodes, total_cells):
        read_count = 0
        read_count_all = 0
        umi_dict = defaultdict(set)
        umi_count = defaultdict()
        with pysam.FastxFile(self.fq2) as fq:
            for read in fq:
                read_count_all+=1
                cb = read.name.split('_')[0]
                umi = read.name.split('_')[1]
                umi_dict[cb].add(umi)
                if cb in cell_barcodes:
                    read_count+=1
        for cb in umi_dict:
            umi_count[cb] = len(umi_dict[cb])
        df_umi = pd.DataFrame.from_dict(umi_count, orient='index', columns=['UMI'])
        df_umi['barcode'] = df_umi.index
        df_umi = df_umi.reset_index(drop=True)
        df_umi = df_umi.reindex(columns=['barcode', 'UMI'])
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in cell_barcodes else 'UB')
        df_umi['barcode'] = df_umi['barcode'].apply(lambda x: self.reversed_compl(x))
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)

        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))
        

        self.summarize_summary.append({
            'item': 'Mean Read Pairs per Cell',
            'count': int(read_count/total_cells),
            'total_count': np.nan
        })

        used_read = 0
        with pysam.FastxFile(self.assembled_fa) as fa:
            for read in fa:
                bc = read.name.split('_')[0]
                if bc in cell_barcodes:
                    used_read += 1
        self.summarize_summary.append({
            'item': 'Mean Used Read Pairs per Cell',
            'count': int(used_read/total_cells), 
            'total_count': np.nan
        })
        self.summarize_summary.append({
            'item': 'Fraction of Reads in Cells',
            'count': used_read,
            'total_count': read_count_all
        })
        
        for c in self.chains:
            temp_df = df_for_clono_pro[df_for_clono_pro['chain']==c]
            self.summarize_summary.append({
                'item': f'Median {c} UMIs per Cell',
                'count': int(temp_df['umis'].median()),
                'total_count': np.nan
            })
        
        annotation_summary = tr.get_vj_annot(df_for_clono, self.chains, self.paired_groups)
        self.summarize_summary += annotation_summary
        # gen stat file          
        stat_file = self.outdir + '/stat.txt'
        sum_df = pd.DataFrame(self.summarize_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)     


    def run(self):
        df_for_clono, df_for_clono_pro, cell_barcodes, total_cells = self.parse_contig_file()
        self.gen_summary(df_for_clono, df_for_clono_pro, cell_barcodes, total_cells)
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