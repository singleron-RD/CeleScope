import pandas as pd
import pysam

from collections import OrderedDict
from celescope.flv_trust4.__init__ import CHAIN, PAIRED_CHAIN
from celescope.flv_CR.summarize import Summarize
from celescope.tools import utils
from celescope.tools.step import s_common, Step
from celescope.tools.plotly_plot import Bar_plot


def gen_vj_annotation_metrics(df, seqtype):
    """
    Generate vdj Annotation Metrics from contig annotations file.
    """
    def get_vj_spanning_pair():
        """
        Get Productive V-J Spanning_Pair metric from annotation file
        Return productive chain pair number. eg: TRA/TRB or IGH/IGL, IGH/IGK.
        """
        df_productive = df[df['productive'] == True]
        
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

    metric_dict = OrderedDict()
    chains, chain_pairs = CHAIN[seqtype], PAIRED_CHAIN[seqtype]

    metric_dict["Cells match with scRNA-seq analysis"] = len(set(df.barcode))
    metric_dict["Cells With Productive V-J Spanning Pair"] = get_vj_spanning_pair()

    for pair in chain_pairs:
        chain1, chain2 = pair.split('_')[0], pair.split('_')[1]
        cbs1 = set(df[(df['productive']==True)&(df['chain']==chain1)].barcode)
        cbs2 = set(df[(df['productive']==True)&(df['chain']==chain2)].barcode)
        metric_dict[f"Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair"] = len(cbs1.intersection(cbs2))

    for chain in chains:
        metric_dict[f"Cells With {chain} Contig"] = len(set(df[df['chain']==chain].barcode))
        metric_dict[f"Cells With CDR3-annotated {chain} Contig"] = len(set(df[(df['chain']==chain)&(df['cdr3']!='None')].barcode))
        metric_dict[f"Cells With V-J Spanning {chain} Contig"] = len(set(df[(df['full_length']==True)&(df['chain']==chain)].barcode))
        metric_dict[f"Cells With Productive {chain} Contig"] =len(set(df[(df['productive']==True)&(df['chain'] == chain)].barcode))

    return metric_dict


def gen_clonotypes_table(df, out_clonotypes):
    """
    Generate clonotypes.csv file
    """
    df = df[df['productive'] == True]
    df['chain_cdr3aa'] = df[['chain', 'cdr3']].apply(':'.join, axis=1)

    match_clonotypes = open(out_clonotypes, 'w')
    match_clonotypes.write('barcode\tcdr3s_aa\n')
    for cb in set(df.barcode):
        temp = df[df['barcode']==cb].sort_values(by='chain', ascending=True)
        chain_pair = ';'.join(temp['chain_cdr3aa'].tolist())
        match_clonotypes.write(f'{cb}\t{chain_pair}\n')
    match_clonotypes.close()

    df_clonetypes = pd.read_csv(out_clonotypes, sep='\t', index_col=None)
    df_clonetypes = df_clonetypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
    df_clonetypes.rename(columns={'barcode': 'frequency'}, inplace=True)
    sum_f = df_clonetypes['frequency'].sum()
    df_clonetypes['proportion'] = df_clonetypes['frequency'].apply(lambda x: x/sum_f)
    df_clonetypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonetypes.shape[0]+1)]
    df_clonetypes = df_clonetypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
    df_clonetypes.sort_values(by='frequency', ascending=False, inplace=True)
    df_clonetypes.to_csv(out_clonotypes, sep=',', index=False)


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
        gen_clonotypes_table(df_match, self.match_clonotypes)

    @utils.add_log
    def gen_matched_metrics(self):
        """
        Generate Matched vdj Annotation Metrics.
        """
        df_match = pd.read_csv(self.match_annotation)
        metrics_dict = gen_vj_annotation_metrics(df_match, self.seqtype)

        for k, v in metrics_dict.items():
            if k == "Cells match with scRNA-seq analysis":
                self.add_metric(
                    name = k,
                    value = v,
                    help_info="The intersection between VDJ cell barcodes and scRNA-Seq barcodes. All the following metrics are based on this intersection."
                )
            else:
                self.add_metric(
                    name = k,
                    value = v,
                    total = len(set(df_match.barcode))
                )

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