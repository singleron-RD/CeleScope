import pandas as pd
import pysam

from celescope.tools import utils
from celescope.tools.step import s_common
from celescope.fl_vdj_CR.VDJ_Mixin import VDJ_Mixin, get_opts_VDJ_Mixin


class Annotation(VDJ_Mixin):
    """
    ## Features

    - V(D)J annotation infomation.

    ## Output
    - `filtered_contig_annotations.csv' High-level annotations of each high-confidence, cellular contig.

    - `filtered_contig.fasta` High-confidence contig sequences in cell barcodes.

    - `clonotypes.csv` High-level descriptions of each clonotype.

    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.seqtype = args.seqtype
        self.barcode_dict = args.barcode_dict

        if self.seqtype == 'TCR':
            self.chains = ['TRA', 'TRB']
            self.pair = ['TRA_TRB']
        elif self.seqtype == 'BCR':
            self.chains = ['IGH', 'IGL', 'IGK']
            self.pair = ['IGH_IGL', 'IGH_IGK']

        self.outs = f'{self.outdir}/../03.assemble/{self.sample}/outs'
        self.all_bam = f'{self.outs}/all_contig.bam'
        self.filter_contig = f'{self.outs}/filtered_contig_annotations.csv'
        self.filter_fa = f'{self.outs}/filtered_contig.fasta'

    @staticmethod
    def get_df_productive(df, seqtype):
        """Get unique productive contig for each chain.

        :param df: filtered contig annotations file.
        :param seqtype: TCR or BCR.
        :return: productive chain pair. eg: TRA/TRB or IGH/IGL, IGH/IGK.
        """
        df_productive = df[(df['full_length'] == True) & (df['productive'] == True)]
        
        if seqtype == "BCR":
            df_chain_heavy = df_productive[ (df_productive['chain'] == 'IGH') ] 
            df_chain_light = df_productive[ (df_productive['chain'] == 'IGL')|(df_productive['chain'] =='IGK')]
            df_chain_heavy.drop_duplicates(['barcode'], inplace=True)
            df_chain_light.drop_duplicates(['barcode'], inplace=True)
            df_productive = pd.merge(df_chain_heavy, df_chain_light, on = 'barcode', how = 'inner')

        else:
            df_TRA = df_productive[ (df_productive['chain'] == 'TRA') ]
            df_TRB = df_productive[ (df_productive['chain'] == 'TRB') ]
            df_TRA.drop_duplicates(['barcode'], inplace=True)
            df_TRB.drop_duplicates(['barcode'], inplace=True)
            df_productive = pd.merge(df_TRA, df_TRB, on = 'barcode', how = 'inner')

        return df_productive
    
    @utils.add_log
    def convert_barcode(self):
        """Convert 10X barcode to SGR barcode

        :return: filter contig annotations in SGR barcode format.
        """

        barcode_df = pd.read_csv(self.barcode_dict, sep='\t', index_col=1)
        barcode_dict = barcode_df.to_dict()['sgr']
        filter_contig = pd.read_csv(f'{self.filter_contig}', sep=',', index_col=None)

        filter_contig['barcode'] = filter_contig['barcode'].apply(lambda x: self.reversed_compl(barcode_dict[x.split('-')[0]]))
        filter_contig['contig_id'] = filter_contig['contig_id'].apply(lambda x: self.reversed_compl(
            barcode_dict[x.split('-')[0]])+'_'+x.split('_')[1]+'_'+x.split('_')[2])
        if self.soft == '3.0.2':
            filter_contig.productive = filter_contig.productive.replace({'True': True, 'None': False})

        filter_contig.to_csv(f'{self.outdir}/filtered_contig_annotations.csv', sep=',', index=False)
    
        in_fa = pysam.FastxFile(self.filter_fa)
        out_fa = open(f'{self.outdir}/filtered_contig.fasta', 'w')
        for entry in in_fa:
            name = entry.name
            seq = entry.sequence
            attrs = name.split('_')
            new_name = self.reversed_compl(barcode_dict[attrs[0].split('-')[0]]) + '_' + attrs[1] + '_' + attrs[2]
            out_fa.write(f'>{new_name}\n{seq}\n')
        out_fa.close()

        return filter_contig
    
    @utils.add_log
    def gen_clonotypes(self):
        """Generate clonotypes.csv file"""

        raw_clonotypes = pd.read_csv(f'{self.outs}/clonotypes.csv', sep=',', index_col=None)
        raw_clonotypes['ClonotypeID'] = raw_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        raw_clonotypes['Frequency'] = raw_clonotypes['frequency']
        raw_clonotypes['Proportion'] = raw_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        raw_clonotypes['CDR3_aa'] = raw_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))
        raw_clonotypes = raw_clonotypes[['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion']]

        raw_clonotypes.to_csv(f'{self.outdir}/clonotypes.csv', sep=',', index=False)

    @utils.add_log
    def get_VDJ_annotation(self, filter_contig, df_productive):
        """VDJ annotation in html"""

        cell_nums = len(set(filter_contig.barcode))

        self.add_metric(
            name = 'Cells With Productive V-J Spanning Pair',
            value = df_productive.shape[0],
            total = cell_nums,
            help_info = "Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
        )

        for pair in self.pair:
            chain1 = pair.split('_')[0]
            chain2 = pair.split('_')[1]
            cbs1 = set(filter_contig[(filter_contig['productive'] == True) & (filter_contig['chain'] == chain1)].barcode)
            cbs2 = set(filter_contig[(filter_contig['productive'] == True) & (filter_contig['chain'] == chain2)].barcode)
            paired_cbs = len(cbs1.intersection(cbs2))

            self.add_metric(
                name = f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
                value = paired_cbs,
                total = cell_nums,
                help_info = "Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
            )

        for chain in self.chains:
            self.add_metric(
                name = f'Cells With {chain} Contig',
                value = len(set(filter_contig[filter_contig['chain'] == chain].barcode)),
                total = cell_nums,
                help_info = f"Fraction of cell-associated barcodes with at least one {chain} contig where a CDR3 was detected"
            )

            self.add_metric(
                name = f'Cells With CDR3-annotated {chain} Contig',
                value = len(set(filter_contig[(filter_contig['chain'] == chain) & (filter_contig['cdr3'] != 'None')].barcode)),
                total = cell_nums,
                help_info = f"Fraction of cell-associated barcodes with at least one {chain} contig where a CDR3 was detected"
            )

            self.add_metric(
                name = f'Cells With V-J Spanning {chain} Contig',
                value = len(set(filter_contig[(filter_contig['full_length'] == True) &
                                 (filter_contig['chain'] == chain)].barcode.tolist())),
                total = cell_nums,
                help_info = f"Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for {chain}"
            )

            self.add_metric(
                name = f'Cells With Productive {chain} Contig',
                value = len(set(filter_contig[(filter_contig['productive'] == True) & (filter_contig['chain'] == chain)].barcode)),
                total = cell_nums,
                help_info = f"raction of cell-associated barcodes with productive {chain} chain. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region"
            )

    def run(self):
        filter_contig = self.convert_barcode()
        self.gen_clonotypes()
        df_productive = self.get_df_productive(filter_contig, self.seqtype)
        self.get_VDJ_annotation(filter_contig, df_productive)


def annotation(args):
    with Annotation(args, display_title="V(D)J Annotation") as runner:
        runner.run()


def get_opts_annotation(parser, sub_program):
    get_opts_VDJ_Mixin(parser)
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--barcode_dict', help='10X barcode correspond sgr barcode', required=True)
    return parser