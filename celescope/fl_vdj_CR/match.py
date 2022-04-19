import pandas as pd
import pysam

from celescope.fl_vdj_CR.annotation import Annotation
from celescope.tools import utils
from celescope.tools.step import s_common
from celescope.fl_vdj_CR.VDJ_Mixin import VDJ_Mixin, get_opts_VDJ_Mixin


class IncorrectMatchDir(Exception):
    pass


class Match(VDJ_Mixin):
    """
    ## Features

    - V(D)J results match SC-RNA infomation.

    ## Output
    - `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

    - `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.
    
    - `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.

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

        self.match_dir = args.match_dir

        self.filter_contig = f'{self.outdir}/../04.annotation/filtered_contig_annotations.csv'
        self.filter_fa = f'{self.outdir}/../03.assemble/{self.sample}/outs/filtered_contig.fasta'

        self.match_fa = f'{self.outdir}/match_contig.fasta'
        
        match_bool = True
        if (not args.match_dir) or (args.match_dir == "None"):
            match_bool = False
        if match_bool:
            try:
                self.match_cell_barcodes, _match_cell_number = utils.get_barcode_from_match_dir(
                    args.match_dir)
            except IndexError as e:
                raise IncorrectMatchDir("Incorrect match_dir, Please Check the match_dir path") from e

    @utils.add_log
    def gen_match_clonotypes(self, df_match):
        """Generate clonotypes file where barcodes match with scRNA

        :param df_match: filter contig annotation where barcodes match with scRNA.
        """
        df_match = df_match[df_match['productive'] == True]
        df_match['chain_cdr3aa'] = df_match[['chain', 'cdr3']].apply(':'.join, axis=1)
        match_cbs = set(df_match.barcode)
        match_clonotypes = open(f'{self.outdir}/match_clonotypes.csv', 'w')
        match_clonotypes.write('barcode\tcdr3s_aa\n')

        for cb in match_cbs:
            temp = df_match[df_match['barcode'] == cb]
            temp = temp.sort_values(by='chain', ascending=True)
            chain_list = temp['chain_cdr3aa'].tolist()
            chain_str = ';'.join(chain_list)
            match_clonotypes.write(f'{cb}\t{chain_str}\n')
        match_clonotypes.close()

        df_match_clonetypes = pd.read_csv(f'{self.outdir}/match_clonotypes.csv', sep='\t', index_col=None)
        df_match_clonetypes = df_match_clonetypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
        df_match_clonetypes = df_match_clonetypes.rename(columns={'barcode': 'frequency'})
        sum_f = df_match_clonetypes['frequency'].sum()
        df_match_clonetypes['proportion'] = df_match_clonetypes['frequency'].apply(lambda x: x/sum_f)
        df_match_clonetypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_match_clonetypes.shape[0]+1)]
        df_match_clonetypes = df_match_clonetypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
        df_match_clonetypes = df_match_clonetypes.sort_values(by='frequency', ascending=False)
        df_match_clonetypes.to_csv(f'{self.outdir}/match_clonotypes.csv', sep=',', index=False)

    @utils.add_log
    def gen_match_fa(self):
        """Generate fasta file where barcodes match with scRNA"""

        barcode_df = pd.read_csv(self.barcode_dict, sep='\t', index_col=1)
        barcode_dict = barcode_df.to_dict()['sgr']
        fa = pysam.FastxFile(self.filter_fa)
        match_fa = open(self.match_fa, 'w')

        for entry in fa:
            name = entry.name
            attrs = name.split('_')
            cb = attrs[0].split('-')[0]
            new_cb = self.reversed_compl(barcode_dict[cb])
            if new_cb in self.match_cell_barcodes:
                new_name = new_cb + '_' + attrs[1] + '_' + attrs[2]
                seq = entry.sequence
                match_fa.write(f'>{new_name}\n{seq}\n')
        match_fa.close() 


    @utils.add_log
    def run(self):
        
        filter_contig = pd.read_csv(self.filter_contig)
        df_RNA_barcode = pd.DataFrame(self.match_cell_barcodes, columns=['barcode'])
        df_match = pd.merge(df_RNA_barcode, filter_contig, on='barcode', how='inner')
        df_match.to_csv(f'{self.outdir}/match_contigs.csv', sep=',', index=False)

        self.add_metric(
            name="Cells match with scRNA-seq analysis",
            value=len(set(df_match.barcode)),
            help_info="Barcodes of cells and barcodes of scRNA-seq cells are reversed complementary"
        )

        df_productive = Annotation.get_df_productive(df_match, self.seqtype)
        Annotation.get_VDJ_annotation(self, df_match, df_productive)
        self.gen_match_clonotypes(df_match)
        self.gen_match_fa()


def match(args):
    with Match(args, display_title="Match") as runner:
        runner.run()


def get_opts_match(parser, sub_program):
    get_opts_VDJ_Mixin(parser)
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--barcode_dict', help='10X barcode correspond sgr barcode', required=True)
        parser.add_argument('--match_dir', help='scRNA-seq match directory', required=True)
    return parser