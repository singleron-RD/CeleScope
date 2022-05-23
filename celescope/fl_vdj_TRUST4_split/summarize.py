from collections import defaultdict
from collections import Counter
import pandas as pd
import pysam
import copy
import os
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.fl_vdj_TRUST4_split.__init__ import CHAIN, PAIRED_CHAIN
from celescope.tools.emptydrop_cr import get_plot_elements


class Summarize(Step):
    """
    ## Features

    - TCR/BCR full length assembly results.

    ## Output
    - `04.summarize/clonetypes.tsv` High-level descriptions of each clonotype.
    - `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
    - `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
    - `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig_annotations.csv.
    - `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
    - `04.summarize/{sample}_two_chain_contig.csv`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.csv.
    - `04.summarize/{sample}_two_chain_contig.fasta`Keep the 2 contigs with the highest UMI. This is a subset of filtered_contig.fasta.
    - `04.summarize/{sample}_one_chain_contig.csv`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.csv.
    - `04.summarize/{sample}_one_chain_contig.fasta`Keep only one chain pair(IGH+IGL/K TRA+TRB) with the highest UMI. This is a subset of chain_filtered_contig.fasta.
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.outdir = args.outdir
        self.sample = args.sample
        self.seqtype = args.seqtype
        self.reads_assignment = args.reads_assignment
        self.fq2 = args.fq2
        self.assembled_fa = args.assembled_fa
        self.record_file = f'{self.outdir}/Cell_num.txt'
        self.coef = int(args.coef)
        self.original_contig = args.contig_file
        self.trust_report = args.trust_report

        self.chains, self.paired_groups = self._parse_seqtype(self.seqtype)

    @staticmethod
    def reversed_compl(seq):
        """Reverse complementary sequence

        :param original seq
        :return Reverse complementary sequence
        """
        return str(Seq(seq).reverse_complement())
    
    @staticmethod
    def _parse_seqtype(seqtype):
        """Parse BCR or TCR

        :param seqtype
        :return CHAIN[seqtype]: 'TCR': ['TRA', 'TRB'], 'BCR': ['IGH', 'IGL', 'IGK']
        	    PAIRED_CHAIN[seqtype]: ['TRA_TRB'], 'BCR': ['IGH_IGL', 'IGH_IGK']
        """
        return CHAIN[seqtype], PAIRED_CHAIN[seqtype]
    
    @staticmethod
    def record_cell_num(df, record_file, step):
        cellnum = len(set(df[df['productive']==True].barcode))
        with open(record_file, 'a+') as f:
            f.write(f'{step}: ' + str(cellnum) + '\n')

    @staticmethod
    def Auto_filter(df, coef):
        """
        Filter contig by assign info of transcriptome, threshold = top 20% contig umi value / coef
        """
        n_cell_1_percentile = len(df) // 5
        sorted_counts = sorted(df.umis, reverse=True)
        count_cell_1_percentile = sorted_counts[n_cell_1_percentile]
        threshold = int(count_cell_1_percentile / coef)
        return threshold

    @utils.add_log
    def parse_contig_file(self):
        """
        Add column name for original contig file.
        """
        df = pd.read_csv(self.original_contig, sep='\t', header=None)
        df.columns = ['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis']
        df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
        df['productive'] = df['cdr3'].apply(lambda x: False if x=='None' else True)

        return df

    @utils.add_log
    def filter_contig(self, df):
        """
        CDR3 filtering
        Filter nonfunctional CDR3(shown 'out_of_frame' in cdr3 report), or CDR3 sequences containing "N" in the nucleotide sequence.
        Keep CDR3aa start with C.
        Keep CDR3aa length >= 5.
        Keep no stop codon in CDR3aa.
        Filter low abundance contigs based on a umi cut-off. 
        """
        df.sort_values(by='umis', ascending=False, inplace=True)
        if self.seqtype == 'BCR':
            df_chain_heavy = df[df['chain']=='IGH']
            df_chain_light = df[(df['chain']=='IGK') | (df['chain']=='IGL')]
        else:
            df_chain_heavy = df[df['chain'] == 'TRA']
            df_chain_light = df[df['chain'] == 'TRB']
        df_chain_heavy = df_chain_heavy.drop_duplicates(['barcode'])
        df_chain_light = df_chain_light.drop_duplicates(['barcode'])
        df_for_clono = pd.concat([df_chain_heavy, df_chain_light], ignore_index=True)
        
        if self.record_file != None:
            Summarize.record_cell_num(df_for_clono, self.record_file, step='Total Assembled Cell Number')
        
        trust_report = pd.read_csv(self.trust_report, sep='\t')
        correct_cdr3 = set(df_for_clono['cdr3']).intersection(set(trust_report.CDR3nt))
        correct_cdr3 = [i for i in correct_cdr3 if i.startswith('C')]
        correct_cdr3 = [i for i in correct_cdr3 if len(i)>=5]
        correct_cdr3 = [i for i in correct_cdr3 if 'UAG' or 'UAA' or 'UGA' not in i]
        df_for_clono = df_for_clono[df_for_clono['cdr3'].isin(correct_cdr3)]

        if self.record_file != None:
            Summarize.record_cell_num(df_for_clono, self.record_file, 
            step='Cell Number: Filter CDR3aa:Not start with C, length<5, no stop codon')
        
        threshold = Summarize.Auto_filter(df_chain_heavy, self.coef)
        df_chain_heavy = df_chain_heavy[df_chain_heavy['umis']>=threshold]
        threshold = Summarize.Auto_filter(df_chain_light, self.coef)
        df_chain_light = df_chain_light[df_chain_light['umis']>=threshold]
        df_for_clono = pd.concat([df_chain_heavy, df_chain_light], ignore_index=True)
        
        if self.record_file != None:
            Summarize.record_cell_num(df_for_clono, self.record_file, 
            step='Cell Number After Auto filter')

        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes = set(df_for_clono_pro['barcode'])

        return df_for_clono, cell_barcodes

    @utils.add_log
    def filter_fasta(self, cell_barcodes):
        """Filter all contig fasta file by barcodes which are identified to be cell.

        :param cell_barcodes: all barcodes identified to be cell.
        """
        all_contig_fasta = f'{self.outdir}/{self.sample}_all_contig.fasta'
        filter_contig_fasta = f'{self.outdir}/{self.sample}_filtered_contig.fasta'

        filter_contig_fasta = open(filter_contig_fasta,'w')
        with pysam.FastxFile(all_contig_fasta) as fa:
            for read in fa:
                name = read.name
                barcode = name.split('_')[0]
                sequence = read.sequence
                if self.reversed_compl(barcode) in cell_barcodes:
                    filter_contig_fasta.write('>' + name + '\n' + sequence + '\n')
                    
        filter_contig_fasta.close()

    @utils.add_log
    def parse_clonotypes(self, df, df_for_clono, cell_barcodes):
        """Parse clonotypes from CDR3 and manually add clonotype id for each contig.

        :param df: original contig file.
        :param df_for_clono: contig info after filter.
        :param cell_barcodes: all barcodes identified to be cell.
        :return df_filter_contig: filtered contigs by cell barcodes.
        """
        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        df_for_clono_pro['chain_cdr3aa'] = df_for_clono_pro[['chain', 'cdr3']].apply(':'.join, axis=1)
        df_for_clono_pro['chain_cdr3nt'] = df_for_clono_pro[['chain', 'cdr3_nt']].apply(':'.join, axis=1)

        cbs = set(df_for_clono_pro['barcode'])
        clonotypes = open(f'{self.outdir}/clonotypes.csv', 'w')
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

        df_clonotypes = pd.read_csv(f'{self.outdir}/clonotypes.csv', sep='\t', index_col=None)
        df_dict = df_clonotypes[["cdr3s_nt", "cdr3s_aa"]].set_index("cdr3s_nt").to_dict(orient='dict')['cdr3s_aa']
        contig_with_clonotype = copy.deepcopy(df_clonotypes)
        df_clonotypes = df_clonotypes.groupby('cdr3s_nt', as_index=False).agg({'barcode': 'count'})
        df_clonotypes.rename(columns={'barcode': 'frequency'}, inplace=True)
        sum_f = df_clonotypes['frequency'].sum()
        df_clonotypes['proportion'] = df_clonotypes['frequency'].apply(lambda x: x/sum_f)
        df_clonotypes.sort_values(by='frequency', ascending=False, inplace=True)
        df_clonotypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonotypes.shape[0]+1)]
        df_clonotypes['cdr3s_aa'] = df_clonotypes['cdr3s_nt'].apply(lambda x:df_dict[x])
        df_clonotypes = df_clonotypes.reindex(columns=['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt'])
        df_clonotypes.to_csv(f'{self.outdir}/clonotypes.csv', sep=',', index=False) 
        used_for_merge = df_clonotypes[['cdr3s_nt','clonotype_id']]

        df_merge = pd.merge(used_for_merge, contig_with_clonotype, on='cdr3s_nt', how='outer')
        df_merge = df_merge[['barcode', 'clonotype_id']]
        df_all_contig = pd.merge(df_merge, df, on='barcode',how='outer')
        df_all_contig.fillna('None',inplace = True)
        df_all_contig = df_all_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_filter_contig = df_all_contig[df_all_contig['barcode'].isin(cell_barcodes)]

        df_all_contig['barcode'] = df_all_contig['barcode'].apply(self.reversed_compl)
        df_all_contig['contig_id'] = df_all_contig['contig_id'].apply(lambda x: self.reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])
        df_all_contig.to_csv(f'{self.outdir}/{self.sample}_all_contig.csv', sep=',', index=False)

        df_filter_contig['barcode'] = df_filter_contig['barcode'].apply(self.reversed_compl)
        df_filter_contig['contig_id'] = df_filter_contig['contig_id'].apply(lambda x: self.reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])
        df_filter_contig.to_csv(f'{self.outdir}/{self.sample}_filtered_contig.csv', sep=',', index=False)

        return df_filter_contig

    @utils.add_log
    def keep_unique_contig(self, df_filter_contig):
        """
        Keep unique contig for each chain(Highest umi) of one cell
        """
        df_filter_contig.sort_values(by='umis',ascending=False)
        if self.seqtype == 'BCR':
            df_chain_heavy = df_filter_contig[df_filter_contig['chain']=='IGH']
            df_chain_light = df_filter_contig[(df_filter_contig['chain']=='IGK') | (df_filter_contig['chain']=='IGL')]
        else:
            df_chain_heavy = df_filter_contig[df_filter_contig['chain'] == 'TRA']
            df_chain_light = df_filter_contig[df_filter_contig['chain'] == 'TRB']
        df_chain_heavy = df_chain_heavy.drop_duplicates(['barcode'])
        df_chain_light = df_chain_light.drop_duplicates(['barcode'])
        df_unique_contig = pd.concat([df_chain_heavy, df_chain_light], ignore_index=True)
        df_unique_contig.to_csv(f'{self.outdir}/{self.sample}_unique_contig.csv', sep=',', index=False)

        unique_contig_set = set(df_unique_contig['contig_id'])
        unique_contig_fasta = f'{self.outdir}/{self.sample}_unique_contig.fasta'
        one_chain_fasta = open(unique_contig_fasta,'w')
        with pysam.FastxFile(f'{self.outdir}/{self.sample}_filtered_contig.fasta') as fa:
            for read in fa:
                seq = read.sequence
                name = read.name
                if name in unique_contig_set:
                    one_chain_fasta.write(">"+name+"\n"+seq+"\n")
        one_chain_fasta.close()

    @utils.add_log
    def gen_summary(self, df_for_clono):
        """ Generate metrics in html 
        """
        df_for_clono_pro = df_for_clono[df_for_clono['productive']==True]
        cell_barcodes = set(df_for_clono_pro['barcode'])
        total_cells =len(cell_barcodes)

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
        df_umi['barcode'] = df_umi['barcode'].apply(self.reversed_compl)
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))

        self.add_metric(
            name = 'Estimated Number of Cells',
            value = total_cells,
            help_info = "Number of cells which contain at least one chain (for TCR: TRA or TRB, for BCR: IGH, IGL or IGK)"            
        )

        self.add_metric(
            name = 'Mean Read Pairs per Cell',
            value = int(read_count/total_cells),
            help_info = 'Number of input reads divided by the estimated number of cells'
        )

        used_read = 0
        with pysam.FastxFile(self.assembled_fa) as fa:
            for read in fa:
                bc = read.name.split('_')[0]
                if bc in cell_barcodes:
                    used_read += 1
        self.add_metric(
            name = 'Mean Used Read Pairs per Cell',
            value = int(used_read/total_cells), 
            help_info = "Mean number of reads used in assembly per cell-associated barcode"
        )
        self.add_metric(
            name = 'Fraction of Reads in Cells',
            value = used_read,
            total = read_count_all,
            help_info = 'Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes'
        )
        
        for c in self.chains:
            temp_df = df_for_clono_pro[df_for_clono_pro['chain']==c]

            try:
                median_umi = int(temp_df['umis'].median())
            except ValueError:
                # ValueError: cannot convert float NaN to integer
                median_umi = 0

            self.add_metric(
                name = f'Median {c} UMIs per Cell',
                value = median_umi
            )

    def run(self):
        original_df = self.parse_contig_file()
        df_for_clono, cell_barcodes = self.filter_contig(original_df)
        self.filter_fasta(cell_barcodes)
        df_filter_contig = self.parse_clonotypes(original_df, df_for_clono, cell_barcodes)
        self.keep_unique_contig(df_filter_contig)
        self.gen_summary(df_for_clono)

@utils.add_log
def summarize(args):
    with Summarize(args, display_title="Cells") as runner:
        runner.run()


def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--coef', help='coef for auto filter', default=10)
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--reads_assignment', help='File records reads assigned to contigs.', required=True)
        parser.add_argument('--fq2', help='Cutadapt R2 reads.', required=True)
        parser.add_argument('--assembled_fa', help='Read used for assembly', required=True)
        parser.add_argument('--trust_report', help='Filtered trust report,Filter Nonfunctional CDR3 and CDR3 sequences containing N', required=True)
        parser.add_argument('--contig_file', help='original contig annotation file', required=True)  
