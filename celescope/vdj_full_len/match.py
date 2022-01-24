import pandas as pd
import pysam
from Bio.Seq import Seq

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.vdj_full_len.annotation import Annotation

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

class Match(Step):
    """
    Features

    - V(D)J results match SC-RNA infomation.

    Output
    - `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

    - `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.
    
    - `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
            # chain parameters
        if self.seqtype == 'TCR':
            self.chains = ['TRA', 'TRB']
            self.pair = ['TRA_TRB']
        elif self.seqtype == 'BCR':
            self.chains = ['IGH', 'IGL', 'IGK']
            self.pair = ['IGH_IGL', 'IGH_IGK']

        self.barcode_dic = args.barcode_dic
        self.match_dir = args.match_dir
        self.remove_3rd_chain = args.remove_3rd_chain

        self.filter_contig = f'{self.outdir}/../04.annotation/filtered_contig_annotations.csv'
        self.filter_fa = f'{self.outdir}/../03.assemble/{self.sample}/outs/filtered_contig.fasta'
        self.match_bool = True
        
        self.match_fa = f'{self.outdir}/match_contig.fasta'
        
        if (not args.match_dir) or (args.match_dir == "None"):
            self.match_bool = False
        if self.match_bool:
            try:
                self.match_cell_barcodes, _match_cell_number = utils.read_barcode_file(
                    args.match_dir)
            except IndexError as e:
                print(
                    "Incorrect match_dir, Please Check the match_dir path" + "\n" + repr(e))
                raise

    @utils.add_log
    def run(self):
        
        filter_contig = pd.read_csv(self.filter_contig)
        df_sgr = pd.DataFrame(self.match_cell_barcodes, columns=['barcode'])
        df_match = pd.merge(df_sgr, filter_contig, on='barcode', how='inner')
        if df_match.empty:
            raise Exception(
                'No match results found in scRNA-seq, please check your match_dir!')

        if self.remove_3rd_chain:
            if self.seqtype == 'BCR':
                df_h = df_match[df_match['chain'] == 'IGH']
                df_temp = df_match[df_match['chain'] != 'IGH']
                df_temp = df_temp.sort_values(by='umis', ascending=False)
                df_temp = df_temp.drop_duplicates(['barcode'])
                df_match = pd.concat([df_h, df_temp], ignore_index=True)

        # get match summary
        self.add_metric(
            name="Cells match with scRNA-seq analysis",
            value=len(set(df_match['barcode'].tolist())),
            help_info="Barcodes of cells and barcodes of scRNA-seq cells are reversed complementary"
        )

        df = df_match
        fl_pro_pair_df = Annotation.Cell_spanning_pair(df, self.seqtype)
        cell_nums = len(set(df['barcode'].tolist()))
        self.add_metric(
            name='Cells With Productive V-J Spanning Pair',
            value=fl_pro_pair_df.shape[0],
            total=cell_nums,
            help_info="Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
        )
        for p in self.pair:
            chain1 = p.split('_')[0]
            chain2 = p.split('_')[1]
            cbs1 = set(df[(df['productive'] == True) & (
                df['chain'] == chain1)].barcode.tolist())
            cbs2 = set(df[(df['productive'] == True) & (
                df['chain'] == chain2)].barcode.tolist())
            paired_cbs = len(cbs1.intersection(cbs2))
            self.add_metric(
                name=f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
                value=paired_cbs,
                total=cell_nums,
                help_info="Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
            )
        for c in self.chains:
            self.add_metric(
                name=f'Cells With {c} Contig',
                value=len(set(df[df['chain'] == c].barcode.tolist())),
                total=cell_nums,
                help_info=f"Fraction of cell-associated barcodes with at least one {c} contig where a CDR3 was detected"
            )
            self.add_metric(
                name=f'Cells With CDR3-annotated {c} Contig',
                value=len(
                    set(df[(df['chain'] == c) & (df['cdr3'] != 'None')].barcode.tolist())),
                total=cell_nums,
                help_info=f"Fraction of cell-associated barcodes with at least one {c} contig where a CDR3 was detected"
            )
            self.add_metric(
                name=f'Cells With V-J Spanning {c} Contig',
                value=len(set(df[(df['full_length'] == True) &
                                 (df['chain'] == c)].barcode.tolist())),
                total=cell_nums,
                help_info=f"Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for {c}"
            )
            self.add_metric(
                name=f'Cells With Productive {c} Contig',
                value=len(set(df[(df['productive'] == True) &
                                 (df['chain'] == c)].barcode.tolist())),
                total=cell_nums,
                help_info=f"raction of cell-associated barcodes with productive {c} chain. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
            )

        df_match.to_csv(
            f'{self.outdir}/match_contigs.csv', sep=',', index=False)

        # get match clonotypes
        df_match = df_match[df_match['productive'] == True]
        df_match['chain_cdr3aa'] = df_match[[
            'chain', 'cdr3']].apply(':'.join, axis=1)
        match_cbs = set(df_match['barcode'].tolist())
        match_clonotypes = open(f'{self.outdir}/match_clonotypes.csv', 'w')
        match_clonotypes.write('barcode\tcdr3s_aa\n')
        for cb in match_cbs:
            temp = df_match[df_match['barcode'] == cb]
            temp = temp.sort_values(by='chain', ascending=True)
            chain_list = temp['chain_cdr3aa'].tolist()
            chain_str = ';'.join(chain_list)
            match_clonotypes.write(f'{cb}\t{chain_str}\n')
        match_clonotypes.close()

        df_match_clonetypes = pd.read_csv(
            f'{self.outdir}/match_clonotypes.csv', sep='\t', index_col=None)
        df_match_clonetypes = df_match_clonetypes.groupby(
            'cdr3s_aa', as_index=False).agg({'barcode': 'count'})
        df_match_clonetypes = df_match_clonetypes.rename(
            columns={'barcode': 'frequency'})
        sum_f = df_match_clonetypes['frequency'].sum()
        df_match_clonetypes['proportion'] = df_match_clonetypes['frequency'].apply(
            lambda x: x/sum_f)
        df_match_clonetypes['clonotype_id'] = [
            f'clonotype{i}' for i in range(1, df_match_clonetypes.shape[0]+1)]
        df_match_clonetypes = df_match_clonetypes.reindex(
            columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
        df_match_clonetypes = df_match_clonetypes.sort_values(
            by='frequency', ascending=False)
        df_match_clonetypes.to_csv(
            f'{self.outdir}/match_clonotypes.csv', sep=',', index=False)

        barcode_df = pd.read_csv(self.barcode_dic, sep='\t', index_col=1)
        barcode_dict = barcode_df.to_dict()['sgr']
        fa = pysam.FastxFile(self.filter_fa)
        match_fa = open(self.match_fa, 'w')
        for entry in fa:
            name = entry.name
            attrs = name.split('_')
            cb = attrs[0].split('-')[0]
            new_cb = reversed_compl(barcode_dict[cb])
            if new_cb in self.match_cell_barcodes:
                new_name = new_cb + '_' + attrs[1] + '_' + attrs[2]
                seq = entry.sequence
                match_fa.write(f'>{new_name}\n{seq}\n')
        match_fa.close()

def match(args):
    with Match(args, display_title="Match") as runner:
        runner.run()

def get_opts_match(parser, sub_program):

    parser.add_argument('--seqtype', help='TCR or BCR',
                        choices=['TCR', 'BCR'], required=True)
    parser.add_argument(
        '--remove_3rd_chain', help='remove IGK or IGL according to umis when a cell has 3 chains at the same time.', action='store_true')
    if sub_program:
        s_common(parser)
        parser.add_argument(
            '--barcode_dic', help='10X barcode correspond sgr barcode', required=True)
        parser.add_argument(
            '--match_dir', help='scRNA-seq match directory', required=True)
    return parser