import pandas as pd
import pysam

from celescope.flv_trust4.__init__ import CHAIN, PAIRED_CHAIN
from celescope.flv_CR.summarize import Summarize
from celescope.tools import utils
from celescope.tools.step import s_common, Step
from celescope.tools.plotly_plot import Bar_plot


class Match(Step):
    """
    ## Features

    - Assembled results match with sc-RNA library.

    - Generate matched VDJ-annotation metrics, clonetypes table and bar-plot of clonetypes distribution in html.

    ## Output

    - `matched_contig_annotations.csv` High-level annotations of each high-confidence contigs from matched cell-associated barcodes.

    - `matched_contig.fasta` High-confidence contig sequences annotated in the matched_contig_annotations.csv.

    - `matched_productive_contig_annotations.csv` Annotations of each productive contigs from matched cell-associated barcodes. This is a subset of matched_contig_annotations.csv.

    - `matched_productive_contig.fasta` Productive contig sequences annotated in the matched_productive_contig_annotations.csv.

    - `clonotypes.csv` High-level descriptions of each clonotype where barcodes match with scRNA-Seq.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.chains = CHAIN[self.seqtype]
        self.pairs = PAIRED_CHAIN[self.seqtype]

        self.match_dir = args.match_dir
        if self.match_dir != 'None':
            self.match_cell_barcodes, _ = utils.get_barcode_from_match_dir(self.match_dir)

        self.filter_annotation = f'{args.summarize_out}/filtered_contig_annotations.csv'
        self.filter_fasta = f'{args.summarize_out}/filtered_contig.fasta'
        self.clonotypes = f'{args.summarize_out}/clonotypes.csv'

        # out
        self.match_annotation = f'{self.outdir}/matched_contig_annotations.csv'
        self.match_fasta = f'{self.outdir}/matched_contig.fasta'
        self.match_clonotypes = f'{self.outdir}/matched_clonotypes.csv'

    @utils.add_log
    def gen_matched_result(self):
        """
        Generate annotation and fasta files where barcodes matched with scRNA.
        """
        SGR_annotation_file = pd.read_csv(self.filter_annotation)
        match_annotation_file = SGR_annotation_file[SGR_annotation_file.barcode.isin(self.match_cell_barcodes)]
        match_annotation_file.to_csv(self.match_annotation, sep=',', index=False)

        SGR_fasta_file = pysam.FastxFile(self.filter_fasta)
        match_fasta_file = open(self.match_fasta, 'w')
        for entry in SGR_fasta_file:
            name = entry.name
            attrs = name.split('_')
            cb = attrs[0]
            if cb in self.match_cell_barcodes:
                new_name = cb + '_' + attrs[1] + '_' + attrs[2]
                seq = entry.sequence
                match_fasta_file.write(f'>{new_name}\n{seq}\n')
        match_fasta_file.close()

        Summarize.gen_productive_contig(match_annotation_file, self.match_fasta, self.outdir, prefix='matched_')
        
    @utils.add_log
    def gen_matched_clonotypes(self):
        """
        Generate clonotypes.csv file where barcodes match with scRNA
        """
        df_match = pd.read_csv(self.match_annotation)
        df_match = df_match[df_match['productive'] == True]
        df_match['chain_cdr3aa'] = df_match[['chain', 'cdr3']].apply(':'.join, axis=1)

        match_clonotypes = open(self.match_clonotypes, 'w')
        match_clonotypes.write('barcode\tcdr3s_aa\n')
        for cb in set(df_match.barcode):
            temp = df_match[df_match['barcode']==cb].sort_values(by='chain', ascending=True)
            chain_pair = ';'.join(temp['chain_cdr3aa'].tolist())
            match_clonotypes.write(f'{cb}\t{chain_pair}\n')
        match_clonotypes.close()

        df_match_clonetypes = pd.read_csv(self.match_clonotypes, sep='\t', index_col=None)
        df_match_clonetypes = df_match_clonetypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
        df_match_clonetypes.rename(columns={'barcode': 'frequency'}, inplace=True)
        sum_f = df_match_clonetypes['frequency'].sum()
        df_match_clonetypes['proportion'] = df_match_clonetypes['frequency'].apply(lambda x: x/sum_f)
        df_match_clonetypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_match_clonetypes.shape[0]+1)]
        df_match_clonetypes = df_match_clonetypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
        df_match_clonetypes.sort_values(by='frequency', ascending=False, inplace=True)
        df_match_clonetypes.to_csv(self.match_clonotypes, sep=',', index=False)

    @utils.add_log
    def gen_matched_metrics(self):
        df_match = pd.read_csv(self.match_annotation)

        self.add_metric(
            name="Cells match with scRNA-seq analysis",
            value=len(set(df_match.barcode)),
            help_info="The intersection between VDJ cell barcodes and scRNA-Seq barcodes. All the following metrics are based on this intersection."
        )

        cell_nums = len(set(df_match.barcode))

        self.add_metric(
            name='Cells With Productive V-J Spanning Pair',
            value=Match.VJ_Spanning_Pair(df_match, self.seqtype),
            total = cell_nums,
        )

        for pair in self.pairs:
            chain1, chain2 = pair.split('_')[0], pair.split('_')[1]
            cbs1 = set(df_match[(df_match['productive']==True)&(df_match['chain']==chain1)].barcode)
            cbs2 = set(df_match[(df_match['productive']==True)&(df_match['chain']==chain2)].barcode)

            self.add_metric(
                name=f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair',
                value=len(cbs1.intersection(cbs2)),
                total=cell_nums,
            )

        for chain in self.chains:
            self.add_metric(
                name=f'Cells With {chain} Contig',
                value=len(set(df_match[df_match['chain']==chain].barcode)),
                total=cell_nums,
            )

            self.add_metric(
                name=f'Cells With CDR3-annotated {chain} Contig',
                value=len(set(df_match[(df_match['chain']==chain)&(df_match['cdr3']!='None')].barcode)),
                total=cell_nums,
            )

            self.add_metric(
                name=f'Cells With V-J Spanning {chain} Contig',
                value=len(set(df_match[(df_match['full_length']==True)&(df_match['chain']==chain)].barcode)),
                total=cell_nums,
            )

            self.add_metric(
                name=f'Cells With Productive {chain} Contig',
                value=len(set(df_match[(df_match['productive']==True)&(df_match['chain'] == chain)].barcode)),
                total=cell_nums,
            )

    @staticmethod
    def VJ_Spanning_Pair(df_annotation, seqtype):
        """
        Get V-J Spanning_Pair metric from annotation file
        Return productive chain pair. eg: TRA/TRB or IGH/IGL, IGH/IGK.
        """
        df_productive = df_annotation[df_annotation['productive'] == True]
        
        if seqtype == "BCR":
            df_chain_heavy = df_productive[(df_productive['chain'] == 'IGH')] 
            df_chain_light = df_productive[(df_productive['chain'] == 'IGL') | (df_productive['chain'] =='IGK')]
        else:
            df_chain_heavy = df_productive[df_productive['chain'] == 'TRA']
            df_chain_light = df_productive[df_productive['chain'] == 'TRB']
        
        for _df in [df_chain_heavy, df_chain_light]:
            _df.drop_duplicates(['barcode'], inplace=True)

        VJ_Spanning_Pair_Cells = pd.merge(df_chain_heavy, df_chain_light, on='barcode', how='inner')
        
        return VJ_Spanning_Pair_Cells.shape[0]

    @utils.add_log
    def gen_clonotypes_table(self):

        title = 'Clonetypes'
        raw_clonotypes = pd.read_csv(self.clonotypes, sep=',', index_col=None)
        raw_clonotypes['ClonotypeID'] = raw_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        raw_clonotypes['Frequency'] = raw_clonotypes['frequency']
        raw_clonotypes['Proportion'] = raw_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        raw_clonotypes['CDR3_aa'] = raw_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

        table_dict = self.get_table_dict(
            title=title,
            table_id='clonetypes',
            df_table=raw_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']]
        )
        self.add_data(table_dict=table_dict)

        raw_clonotypes['ClonotypeID'] = raw_clonotypes['ClonotypeID'].astype("int")
        raw_clonotypes.sort_values(by=['ClonotypeID'], inplace=True)
        Barplot = Bar_plot(df_bar=raw_clonotypes).get_plotly_div()
        self.add_data(Barplot=Barplot)

    @utils.add_log
    def run(self):
        if self.match_dir != 'None':
            self.gen_matched_result()
            self.gen_matched_clonotypes()
            self.gen_matched_metrics()
        self.gen_clonotypes_table()


def match(args):
    with Match(args, display_title="Match") as runner:
        runner.run()


def get_opts_match(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
        parser.add_argument('--summarize_out', help='assemble result in SGR barcode from summarize directory', required=True)
    return parser