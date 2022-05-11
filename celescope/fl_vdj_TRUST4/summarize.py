from collections import defaultdict
import pandas as pd
import pysam
import copy
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.fl_vdj_TRUST4.__init__ import CHAIN, PAIRED_CHAIN
from celescope.tools.emptydrop_cr import get_plot_elements


@utils.add_log
def gen_contig_csv(outdir, sample, full_len_assembly, assign_out):
    """ Generate contig file in 10X format.
    
    Genrate contig file from full length annotated fasta file.

    Args:
        outdir: output directory.
        sample: sample name.
        full_len_assembly: full length annotated fasta file.
        assign_out: read assignment results.
    """
    # reads assignment 
    assignment = pd.read_csv(assign_out, sep='\t', header=None)
    # assignment['read_barcode'] = assignment[0].apply(lambda x: x.split('_')[0])
    # assignment['contig_barcode'] = assignment[1].apply(lambda x: x.split('_')[0])
    # assignment['match_barcode'] = assignment[['read_barcode', 'contig_barcode']].apply(lambda x: x['read_barcode']==x['contig_barcode'], axis=1)
    # assignment = assignment[assignment['match_barcode']==True]
    # assignment['umi'] = assignment[0].apply(lambda x: x.split('_')[1])

    assignment = assignment.rename(columns={0:'read_name',1:'contig_id'})
    assignment['umi'] = assignment['read_name'].apply(lambda x:x.split('_')[1])
    read_count_dict = assignment.groupby('contig_id')['read_name'].apply(lambda x:len(set(x))).to_dict()
    umi_count_dict = assignment.groupby('contig_id')['umi'].apply(lambda x:len(set(x))).to_dict()
    # write contig csv
    contigs = open(f'{outdir}/{sample}_contig.csv', 'w')
    process_read = 0
    with pysam.FastxFile(full_len_assembly) as fa:
        for read in fa:
            name = read.name
            comment = read.comment
            attrs = comment.split(' ')
            barcode = name.split('_')[0]
            is_cell = 'True'
            high_confidence = 'True'
            length = attrs[0]
            chain = attrs[2][:3]
            full_length = 'True'
            v_gene = attrs[2].split('(')[0]
            d_gene = attrs[3]
            j_gene = attrs[4].split('(')[0]
            c_gene = attrs[5]
            cdr3 = attrs[8].split('=')[1]
            cdr3_aa = 'None'
            productive = 'False'
            # temp = assignment[assignment[1]==name]
            # reads = str(len(temp[0].tolist()))
            # umis = str(len(set(temp['umi'].tolist())))
            reads = str(read_count_dict.get(name, 0))
            umis = str(umi_count_dict.get(name, 0))

            string = '\t'.join([barcode, is_cell, name, high_confidence, length, chain, v_gene, d_gene, j_gene, c_gene, full_length, productive, cdr3_aa, cdr3, reads, umis])
            contigs.write(f'{string}\n')
            process_read+=1
            if process_read % 10000 == 0:
                gen_contig_csv.logger.info(f'Processed {process_read} contigs')

    contigs.close()


class Summarize(Step):

    """
    ## Features

    - TCR/BCR full length assembly results.

    ## Output
    - `04.summarize/clonetypes.csv` High-level descriptions of each clonotype.
    - `04.summarize/{sample}_all_contig.csv` High-level and detailed annotations of each contig.
    - `04.summarize/{sample}_all_contig.fasta` All assembled contig sequences.
    - `04.summarize/{sample}_filtered_contig.csv` High-level annotations of each cellular contig after filter. This is a subset of all_contig.csv.
    - `04.summarize/{sample}_filtered_contig.fasta` Assembled contig sequences after filter.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.full_len_assembly = args.full_len_assembly
        self.assign_out = args.assign_out
        self.filter_report = args.filter_report
        self.barcode_filter_report = args.barcode_filter_report
        self.cutadapted_fq = args.cutadapted_fq
        self.assembled_fa = args.assembled_fa

        self.chains, self.paired_groups = self._parse_seqtype(self.seqtype)

    @staticmethod
    def reversed_compl(seq):
        """ Return reverse complementation sequence."""
        return str(Seq(seq).reverse_complement())
    
    @staticmethod
    def _parse_seqtype(seqtype):
        """ Return 'TRA_TRB' for TCR or 'IGH_IGL', 'IGH_IGK' for BCR."""
        return CHAIN[seqtype], PAIRED_CHAIN[seqtype]
    
    @staticmethod
    def cut_off(barcode_count):
        """ barcode cut-off by umi and reads.
        
        Args:
            barcode_count: umi and read count info for each barcode.
        
        Returns:
            barcode_count: filtered barcode_count dataframe.
        """
        barcode_count.sort_values('umis', ascending=True, inplace = True)
        RANK = barcode_count.shape[0]//5
        rank_UMI = barcode_count.iloc[RANK, :]["umis"]
        UMI_min = int(rank_UMI)
        
        barcode_count.sort_values('reads', ascending=True, inplace = True)
        RANK = barcode_count.shape[0]//5
        rank_read = barcode_count.iloc[RANK, :]["reads"]
        UMI_read = int(rank_read)

        barcode_count = barcode_count[barcode_count['umis']>=UMI_min]
        barcode_count = barcode_count[barcode_count['reads']>=UMI_read]

        return barcode_count

    @utils.add_log
    def gen_all_contig_file(self):
        """ Generate sequence info of all assembled contigs in fasta format."""
        all_fa = open(f'{self.outdir}/{self.sample}_all_contig.fasta','w')
        with pysam.FastxFile(self.full_len_assembly) as fa:
            for read in fa: 
                name = read.name
                barcode = name.split('_')[0]
                sequence = read.sequence
                all_fa.write('>' + self.reversed_compl(barcode) + '_' + name.split('_')[1] + '\n' + sequence + '\n')    
        all_fa.close()

        gen_contig_csv(self.outdir, self.sample, self.full_len_assembly, self.assign_out)
    
    @utils.add_log
    def parse_contig_file(self):
        """ Return formatted contig file."""
        df = pd.read_csv(f'{self.outdir}/{self.sample}_contig.csv', sep='\t', header=None)
        df.columns = ['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis']
        df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
        df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
        df['productive'] = df['cdr3'].apply(lambda x: False if x=='None' else True)

        return df

    @staticmethod
    def filter_cell(df, seqtype, filter_report, barcode_filter_report):
        """ Filter barcode in contig file.
        Keep CDR3aa start with C.
        Keep CDR3aa length >= 5.
        Keep no stop codon in CDR3aa.
        Filter low abundance contigs in barcode.

        Args:
            df: contig file.
            seqtype: BCR or TCR.
            filter_report: filtered trust report file.
            barcode_filter_report: filtered barcode report file.
        
        Returns:
            df_filter: filtered contig info.
            productive_barcodes: productive barcode set in filtered contig file.
        """
        filter_report = pd.read_csv(filter_report, sep='\t')
        barcode_filter_report = pd.read_csv(barcode_filter_report, sep='\t')

        df.sort_values(by='umis', ascending=False, inplace=True)
        if seqtype == 'BCR':
            df_chain_heavy = df[df['chain']=='IGH']
            df_chain_light = df[(df['chain']=='IGK') | (df['chain']=='IGL')]
            df_chain_heavy = df_chain_heavy.drop_duplicates(['barcode'])
            df_chain_light = df_chain_light.drop_duplicates(['barcode'])
            df_filter = pd.concat([df_chain_heavy, df_chain_light], ignore_index=True)
        else:
            df_TRA = df[df['chain'] == 'TRA']
            df_TRB = df[df['chain'] == 'TRB']
            df_TRA = df_TRA.drop_duplicates(['barcode'])
            df_TRB = df_TRB.drop_duplicates(['barcode'])
            df_filter = pd.concat([df_TRA, df_TRB], ignore_index=True)
        
        barcode_filter_report.rename(columns = {'#barcode':'barcode'}, inplace=True)

        filter_report = filter_report[filter_report['cid_full_length'] >= 1]
        filter_report.rename(columns = {'cid':'barcode', '#count':'count'}, inplace=True)
        filter_report['barcode'] = filter_report['barcode'].apply(lambda x:x.split('_')[0])
        cdr3_list = list(filter_report['CDR3aa'])
        cdr3_list = [i for i in cdr3_list if i.startswith('C')]
        cdr3_list = [i for i in cdr3_list if len(i)>=5]
        cdr3_list = [i for i in cdr3_list if 'UAG' or 'UAA' or 'UGA' not in i]
        filter_report = filter_report[filter_report['CDR3aa'].isin(cdr3_list)]

        # df_filter = df_filter[df_filter['umis']>=3]
        # df_filter = df_filter[df_filter['reads']>=2]
        df_filter = df_filter[df_filter['barcode'].isin(set(filter_report['barcode']))]
        df_filter = df_filter[df_filter['barcode'].isin(set(barcode_filter_report['barcode']))]

        # record file IGH+IGK/IGL for BCR, TRA+TRB for TCR
        barcode_count = df_filter.groupby(['barcode']).agg({'umis': 'mean','reads': 'mean'}).reset_index()
        filter_barcode_count = Summarize.cut_off(barcode_count)

        df_filter = df_filter[df_filter['barcode'].isin(filter_barcode_count.barcode)]

        # df_chain_pair.to_csv(f'{self.outdir}/{self.sample}_chain_pair.csv', sep=',', index=False)
        df_filter_pro = df_filter[df_filter['productive']==True]
        productive_barcodes = set(df_filter_pro['barcode'])

        return df_filter, productive_barcodes

    @utils.add_log
    def gen_filter_fasta(self, productive_barcodes):
        """ Generate filtered contig fasta file."""
        all_contig_fa = f'{self.outdir}/{self.sample}_all_contig.fasta'
        out_filter_fa = open(f'{self.outdir}/{self.sample}_filtered_contig.fasta','w')
        with pysam.FastxFile(all_contig_fa) as fa:
            for read in fa:
                name = read.name
                barcode = name.split('_')[0]
                sequence = read.sequence
                if Summarize.reversed_compl(barcode) in productive_barcodes:
                    out_filter_fa.write('>' + name + '\n' + sequence + '\n')
        out_filter_fa.close()
    
    @staticmethod
    def reverse_contig(df, productive_barcodes):
        """ Return dataframe of all contig and filtered contig that including reverse complementation barcode info."""
        # all contig.csv
        df_all_contig = copy.deepcopy(df)
        df_all_contig['barcode'] = df_all_contig['barcode'].apply(Summarize.reversed_compl)
        df_all_contig['contig_id'] = df_all_contig['contig_id'].apply(lambda x: Summarize.reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])
        # filter contig.csv
        df_filter_contig = copy.deepcopy(df)
        df_filter_contig = df_filter_contig[df_filter_contig['barcode'].isin(productive_barcodes)]
        df_filter_contig['barcode'] = df_filter_contig['barcode'].apply(Summarize.reversed_compl)
        df_filter_contig['contig_id'] = df_filter_contig['contig_id'].apply(lambda x: Summarize.reversed_compl(x.split('_')[0]) + '_' + x.split('_')[1])

        return df_all_contig, df_filter_contig

    @utils.add_log
    def gen_clonotypes_file(self, df_chain_pair, productive_barcodes):
        """ Generate clonotypes.csv file.

        Args:
            df_chain_pair: filtered contig info.
            productive_barcodes: productive barcode set in filtered contig file.
        Returns:
            contig_with_clonotype: dataframe records barcode, cdr3s_aa, cdr3s_nt info.
            used_for_merge: datafram records cdr3s_nt, clonotype_id info.
        """
        df_chain_pair['chain_cdr3aa'] = df_chain_pair[['chain', 'cdr3']].apply(':'.join, axis=1)
        df_chain_pair['chain_cdr3nt'] = df_chain_pair[['chain', 'cdr3_nt']].apply(':'.join, axis=1)

        clonotypes = open(f'{self.outdir}/clonotypes.csv', 'w')
        clonotypes.write('barcode\tcdr3s_aa\tcdr3s_nt\n')
        for cb in productive_barcodes:
            temp = df_chain_pair[df_chain_pair['barcode']==cb]
            temp = temp.sort_values(by='chain', ascending=True)
            chain_list = temp['chain_cdr3aa'].tolist()
            chain_str = ';'.join(chain_list)
            ntchain_list = temp['chain_cdr3nt'].tolist()
            ntchain_str = ';'.join(ntchain_list)
            clonotypes.write(f'{cb}\t{chain_str}\t{ntchain_str}\n')
        clonotypes.close()

        df_clonotypes = pd.read_csv(f'{self.outdir}/clonotypes.csv', sep='\t', index_col=None)
        nt_aa_dict = df_clonotypes[["cdr3s_nt", "cdr3s_aa"]].set_index("cdr3s_nt").to_dict(orient='dict')['cdr3s_aa']
        contig_with_clonotype = copy.deepcopy(df_clonotypes)

        df_clonotypes = df_clonotypes.groupby('cdr3s_nt', as_index=False).agg({'barcode': 'count'})
        df_clonotypes = df_clonotypes.rename(columns={'barcode': 'frequency'})
        sum_frequency = df_clonotypes['frequency'].sum()
        df_clonotypes['proportion'] = df_clonotypes['frequency'].apply(lambda x: x/sum_frequency)
        df_clonotypes = df_clonotypes.sort_values(by='frequency', ascending=False)
        df_clonotypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonotypes.shape[0]+1)]
        df_clonotypes['cdr3s_aa'] = df_clonotypes['cdr3s_nt'].apply(lambda x:nt_aa_dict[x])
        df_clonotypes = df_clonotypes.reindex(columns=['clonotype_id', 'frequency', 'proportion', 'cdr3s_aa', 'cdr3s_nt'])
        df_clonotypes.to_csv(f'{self.outdir}/clonotypes.csv', sep=',', index=False) 
        used_for_merge = df_clonotypes[['cdr3s_nt','clonotype_id']]

        return contig_with_clonotype, used_for_merge

    @utils.add_log
    def add_clonotypes_for_contig(self, used_for_merge, contig_with_clonotype, df_all_contig, df_filter_contig):
        """ Add clonotype_id in contig.csv file.
        
        Args:
            contig_with_clonotype: dataframe records barcode, cdr3s_aa, cdr3s_nt info.
            used_for_merge: datafram records cdr3s_nt, clonotype_id info.
            df_all_contig: all contig info.
            df_filter_contig: filtered contig info.
        """
        df_merge = pd.merge(used_for_merge, contig_with_clonotype, on='cdr3s_nt', how='outer')
        df_merge = df_merge[['barcode','clonotype_id']]
        df_merge['barcode'] = df_merge['barcode'].apply(Summarize.reversed_compl)
        df_all_contig = pd.merge(df_merge, df_all_contig, on='barcode',how='outer')
        df_filter_contig = pd.merge(df_merge, df_filter_contig, on='barcode',how='outer')
        df_all_contig.fillna('None',inplace = True)
        df_filter_contig.fillna('None',inplace = True)
        df_all_contig = df_all_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_filter_contig = df_filter_contig[['barcode', 'is_cell', 'contig_id', 'high_confidence', 'length', 'chain', 'v_gene', 'd_gene', 'j_gene', 'c_gene', 'full_length', 'productive', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'clonotype_id']]
        df_all_contig.to_csv(f'{self.outdir}/{self.sample}_all_contig.csv', sep=',', index=False)
        df_filter_contig.to_csv(f'{self.outdir}/{self.sample}_filtered_contig.csv', sep=',', index=False)
    
    @utils.add_log
    def gen_summary(self, df_chain_pair, productive_barcodes):
        """ Generate summary metrics for html.
        Including: Estimated Number of Cells, Mean Read Pairs per Cell, Mean Used Read Pairs per Cell,
            Fraction of Reads in Cells, Median UMIs per Cell, Barcode rank plot

        Args:
            df_chain_pair: filtered contig info
            productive_barcodes: productive barcode set in filtered contig.
        """
        total_cell_num = len(productive_barcodes)
        read_count = 0
        read_count_all = 0
        umi_dict = defaultdict(set)
        umi_count = defaultdict()
        with pysam.FastxFile(self.cutadapted_fq) as fq:
            for read in fq:
                read_count_all+=1
                cb = read.name.split('_')[0]
                umi = read.name.split('_')[1]
                umi_dict[cb].add(umi)
                if cb in productive_barcodes:
                    read_count+=1
        for cb in umi_dict:
            umi_count[cb] = len(umi_dict[cb])
        df_umi = pd.DataFrame.from_dict(umi_count, orient='index', columns=['UMI'])
        df_umi['barcode'] = df_umi.index
        df_umi = df_umi.reset_index(drop=True)
        df_umi = df_umi.reindex(columns=['barcode', 'UMI'])
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in productive_barcodes else 'UB')
        df_umi['barcode'] = df_umi['barcode'].apply(self.reversed_compl)
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))

        self.add_metric(
            name = 'Estimated Number of Cells',
            value = total_cell_num,
            help_info = "Number of cells which contain at least one chain (for TCR: TRA or TRB, for BCR: IGH, IGL or IGK)"            
        )

        self.add_metric(
            name = 'Mean Read Pairs per Cell',
            value = int(read_count/total_cell_num),
            help_info = 'Number of input reads divided by the estimated number of cells'
        )

        used_read = 0
        with pysam.FastxFile(self.assembled_fa) as fa:
            for read in fa:
                bc = read.name.split('_')[0]
                if bc in productive_barcodes:
                    used_read += 1
        self.add_metric(
            name = 'Mean Used Read Pairs per Cell',
            value = int(used_read/total_cell_num), 
            help_info = "Mean number of reads used in assembly per cell-associated barcode"
        )
        self.add_metric(
            name = 'Fraction of Reads in Cells',
            value = used_read,
            total = read_count_all,
            help_info = 'Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes'
        )
        
        for c in self.chains:
            temp_df = df_chain_pair[df_chain_pair['chain']==c]
            self.add_metric(
                name = f'Median {c} UMIs per Cell',
                value = int(temp_df['umis'].median())
            )

    def run(self):
        self.gen_all_contig_file()
        df = self.parse_contig_file()
        df_chain_pair, productive_barcodes = self.filter_cell(df, self.seqtype, self.filter_report, self.barcode_filter_report)
        self.gen_filter_fasta(productive_barcodes)
        df_all_contig, df_filter_contig = self.reverse_contig(df, productive_barcodes)
        contig_with_clonotype, used_for_merge = self.gen_clonotypes_file(df_chain_pair, productive_barcodes)
        self.add_clonotypes_for_contig(used_for_merge, contig_with_clonotype, df_all_contig, df_filter_contig)
        self.gen_summary(df_chain_pair, productive_barcodes)


@utils.add_log
def summarize(args):
    with Summarize(args, display_title="Cells") as runner:
        runner.run()


def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--full_len_assembly', help='full len assembly file', required=True)
        parser.add_argument('--assign_out', help='read assignment results', required=True)
        parser.add_argument('--filter_report', help='Filtered trust report', required=True)
        parser.add_argument('--barcode_filter_report', help='Filtered barcode report', required=True)
        parser.add_argument('--cutadapted_fq', help='Cutadapt R2 reads.', required=True)
        parser.add_argument('--assembled_fa', help='Read used for assembly', required=True)
