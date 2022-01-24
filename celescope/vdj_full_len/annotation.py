import pandas as pd
import pysam
from Bio.Seq import Seq

from celescope.tools import utils
from celescope.tools.step import Step, s_common

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

class Annotation(Step):
    """
    Features

    - V(D)J annotation infomation.

    Output
    - `filtered_contig_annotations.csv' 

    - `filtered_contig.fasta`

    - `clonotypes.csv` 

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.soft = args.soft
        self.seqtype = args.seqtype
            # chain parameters
        if self.seqtype == 'TCR':
            self.chains = ['TRA', 'TRB']
            self.pair = ['TRA_TRB']
        elif self.seqtype == 'BCR':
            self.chains = ['IGH', 'IGL', 'IGK']
            self.pair = ['IGH_IGL', 'IGH_IGK']

        # input
        self.outs = f'{self.outdir}/../03.assemble/{self.sample}/outs'
        self.all_bam = f'{self.outs}/all_contig.bam'
        self.filter_contig = f'{self.outs}/filtered_contig_annotations.csv'
        self.filter_fa = f'{self.outs}/filtered_contig.fasta'
        self.barcode_dic = args.barcode_dic

    @staticmethod
    def Cell_spanning_pair(df, seqtype):
        if seqtype == "BCR":
            fl_pro_pair_df = df[(df['full_length'] == True) & (df['productive'] == True)]
            df_chain_heavy = fl_pro_pair_df[ (fl_pro_pair_df['chain'] == 'IGH') ] 
            df_chain_light = fl_pro_pair_df[ (fl_pro_pair_df['chain'] == 'IGL')|(fl_pro_pair_df['chain'] =='IGK')]
            df_chain_heavy = df_chain_heavy.drop_duplicates(['barcode'])
            df_chain_light = df_chain_light.drop_duplicates(['barcode'])
            fl_pro_pair_df = pd.merge(df_chain_heavy, df_chain_light, on = 'barcode', how = 'inner')
        else:
            fl_pro_pair_df = pd.DataFrame(df[(df['full_length'] == True) & (
                df['productive'] == True)].drop_duplicates(['barcode', 'chain']).barcode.value_counts())
            fl_pro_pair_df = fl_pro_pair_df[fl_pro_pair_df['barcode'] >= 2]
        return fl_pro_pair_df
    
    @utils.add_log
    def run(self):
        barcode_df = pd.read_csv(self.barcode_dic, sep='\t', index_col=1)
        barcode_dict = barcode_df.to_dict()['sgr']

        # convert barcode to sgr_cbs, filtered_contig_annotations.csv
        filter_contig = pd.read_csv(
            f'{self.filter_contig}', sep=',', index_col=None)
        filter_contig['barcode'] = filter_contig['barcode'].apply(
            lambda x: reversed_compl(barcode_dict[x.split('-')[0]]))
        filter_contig['contig_id'] = filter_contig['contig_id'].apply(lambda x: reversed_compl(
            barcode_dict[x.split('-')[0]])+'_'+x.split('_')[1]+'_'+x.split('_')[2])
        if self.soft == '3.0.2':
            filter_contig.productive = filter_contig.productive.replace(
                {'True': True, 'None': False})
        filter_contig.to_csv(
            f'{self.outdir}/filtered_contig_annotations.csv', sep=',', index=False)
        
        # filtered_contig.fasta
        in_fa = pysam.FastxFile(self.filter_fa)
        out_fa = open(f'{self.outdir}/filtered_contig.fasta', 'w')
        for entry in in_fa:
            name = entry.name
            seq = entry.sequence
            attrs = name.split('_')
            new_name = reversed_compl(barcode_dict[attrs[0].split(
                '-')[0]]) + '_' + attrs[1] + '_' + attrs[2]
            out_fa.write(f'>{new_name}\n{seq}\n')
        out_fa.close()

        # gen clone table
        raw_clonotypes = pd.read_csv(
            f'{self.outs}/clonotypes.csv', sep=',', index_col=None)
        raw_clonotypes['ClonotypeID'] = raw_clonotypes['clonotype_id'].apply(
            lambda x: x.strip('clonetype'))
        raw_clonotypes['Frequency'] = raw_clonotypes['frequency']
        raw_clonotypes['Proportion'] = raw_clonotypes['proportion'].apply(
            lambda x: f'{round(x*100, 2)}%')
        raw_clonotypes['CDR3_aa'] = raw_clonotypes['cdr3s_aa'].apply(
            lambda x: x.replace(';', '<br>'))
        all_clonotypes = raw_clonotypes[[
            'clonotype_id', 'cdr3s_aa', 'frequency', 'proportion']]
        all_clonotypes.to_csv(
            f'{self.outdir}/clonotypes.csv', sep=',', index=False)

        # V(D)J Annotation summary
        df = filter_contig
        fl_pro_pair_df = self.Cell_spanning_pair(df, self.seqtype)
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

def annotation(args):
    with Annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()

def get_opts_annotation(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR',
                        choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--soft', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'],
                        default='4.0.0')
    if sub_program:
        s_common(parser)
        parser.add_argument(
            '--barcode_dic', help='10X barcode correspond sgr barcode', required=True)
    return parser