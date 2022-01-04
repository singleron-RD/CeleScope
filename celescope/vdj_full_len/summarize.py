import os
from collections import defaultdict
import pandas as pd
import pysam
from Bio.Seq import Seq

from celescope.tools import utils
from celescope.tools import step
from celescope.tools.cellranger3 import get_plot_elements
from celescope.tools.step import Step, s_common


def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())


class Summarize(Step):
    """
    Features

    - Summarize contig and clonetypes infomation.

    Output
    - `filtered_contig_annotations.csv' annotations of each high-confidence, cellular contig.

    - `filtered_contig.fasta` The clone type of VDJ cell barcode, each cell occupies one line.

    - `clonotypes.csv` The count and percentage of each clonetypes of VDJ cell barcode.

    - `match_contigs.csv` filtered_contig_annotations.csv matched with scRNA-seq.

    - `match_contig.fasta` filtered_contig.fasta matched with scRNA-seq.
    
    - `match_clonotypes.csv` clonetypes result matched with scRNA-seq

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.barcode_dic = args.barcode_dic
        self.all_bam = f'{self.outdir}/../03.assemble/{self.sample}/outs/all_contig.bam'
        self.annotation = f'{self.outdir}/../04.annotation'
        self.match = f'{self.outdir}/../05.match'
        self.filter_contig = pd.read_csv(f'{self.outdir}/../04.annotation/filtered_contig_annotations.csv')
        self.clonotable = f'{self.outdir}/../03.assemble/{self.sample}/outs/clonotypes.csv'

        if self.seqtype == 'TCR':
            self.chains = ['TRA', 'TRB']
            self.pair = ['TRA_TRB']
        elif self.seqtype == 'BCR':
            self.chains = ['IGH', 'IGL', 'IGK']
            self.pair = ['IGH_IGL', 'IGH_IGK']


    @utils.add_log
    def run(self):

        # copy file
        os.system(f'cp {self.annotation}/* {self.outdir}')
        os.system(f'cp {self.match}/* {self.outdir}')

        # get barcode correspond dictionary
        barcode_df = pd.read_csv(self.barcode_dic, sep='\t', index_col=1)
        barcode_dict = barcode_df.to_dict()['sgr']

        # get barcode_stat
        stat_dict = pd.read_csv(
            f'{self.outdir}/../01.barcode/stat.txt', sep=':', header=None)
        read_count = int(stat_dict.iloc[0, 1].replace(',', ''))
        sum_dict = pd.read_csv(
            f'{self.outdir}/../03.assemble/{self.sample}/outs/metrics_summary.csv', sep=',', index_col=None)
        sum_dict = sum_dict.T.to_dict()
        total_reads = int(sum_dict[0]["Number of Read Pairs"].replace(',', ''))


        if self.seqtype == "TCR":
            cell_nums = len(set(self.filter_contig['barcode'].tolist()))
            self.add_metric(
                name='Estimated Number of Cells',
                value=cell_nums,
                help_info="Cells with at least one productive TRA or TRB chain"
            )
            self.add_metric(
                name='Reads Mapped to Any V(D)J Gene',
                value=int(
                    total_reads * (float(sum_dict[0]['Reads Mapped to Any V(D)J Gene'].strip('%'))/200)),
                total=read_count,
                help_info="Reads mapped to any TCR or BCR genes."
            )
            for c in self.chains:
                self.add_metric(
                    name=f'Reads Mapped to {c}',
                    value=int(
                        total_reads * (float(sum_dict[0][f'Reads Mapped to {c}'].strip('%'))/200)),
                    total=read_count,
                    help_info=f"Reads mapped to {c} chain. For BCR, this should be one of [IGH, IGL, IGK]"
                )

            self.add_metric(
                name='Fraction Reads in Cells',
                value=int(
                    total_reads * (float(sum_dict[0]['Fraction Reads in Cells'].strip('%'))/200)),
                total=read_count,
                help_info="Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes"
            )
            for c in self.chains:
                mid = self.filter_contig[self.filter_contig['chain']
                                    == c]['umis'].median()
                if mid == mid:
                    self.add_metric(
                        name=f'Median used {c} UMIs per Cell',
                        value=int(mid),
                        help_info=f"Median number of UMIs assigned to a {c} contig per cell. For B cells, only the max of [IGK, IGL] are counted"
                    )
                else:
                    self.add_metric(
                        name=f'Median used {c} UMIs per Cell',
                        value=0,
                        help_info=f"Median number of UMIs assigned to a {c} contig per cell. For B cells, only the max of [IGK, IGL] are counted"
                    )

        elif self.seqtype == 'BCR':
            cell_nums = len(set(self.filter_contig['barcode'].tolist()))
            self.add_metric(
                name='Estimated Number of Cells',
                value=cell_nums,
                help_info="Cells with at least one productive IGH, IGL or IGK chain"
            )
            self.add_metric(
                name='Reads Mapped to Any V(D)J Gene',
                value=int(
                    total_reads * (float(sum_dict[0]['Reads Mapped to Any V(D)J Gene'].strip('%'))/200)),
                total=read_count,
                help_info="Reads mapped to any IGH, IGL or IGK genes"
            )
            for c in self.chains:
                self.add_metric(
                    name=f'Reads Mapped to {c}',
                    value=int(
                        total_reads * (float(sum_dict[0][f'Reads Mapped to {c}'].strip('%'))/200)),
                    total=read_count,
                    help_info=f"Reads mapped to {c} chain. For BCR, this should be one of [TRA, TRB]"
                )

            self.add_metric(
                name='Fraction Reads in Cells',
                value=int(
                    total_reads * (float(sum_dict[0]['Fraction Reads in Cells'].strip('%'))/200)),
                total=read_count,
                help_info="Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes"
            )
            for c in self.chains:
                mid = self.filter_contig[self.filter_contig['chain']
                                    == c]['umis'].median()
                if mid == mid:
                    self.add_metric(
                        name=f'Median used {c} UMIs per Cell',
                        value=int(mid),
                        help_info=f"Median number of UMIs assigned to a {c} contig per cell."
                    )
                else:
                    self.add_metric(
                        name=f'Median used {c} UMIs per Cell',
                        value=0,
                        help_info=f"Median number of UMIs assigned to a {c} contig per cell."
                    )

        # count umi and plot
        all_bam = pysam.AlignmentFile(self.all_bam)
        dic_umi = defaultdict(set)

        for read in all_bam:
            cb = read.get_tag('CB')
            umi = read.get_tag('UB')
            new_cb = reversed_compl(barcode_dict[cb.split('-')[0]])
            dic_umi[new_cb].add(umi)

        df_umi = pd.DataFrame()
        df_umi['barcode'] = list(dic_umi.keys())
        df_umi['UMI'] = [len(dic_umi[i]) for i in dic_umi]
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        sgr_cbs = set(self.filter_contig['barcode'].tolist())
        df_umi['mark'] = df_umi['barcode'].apply(
            lambda x: 'CB' if x in sgr_cbs else 'UB')
        # self.add_data(CB_num=df_umi[df_umi['mark'] == 'CB'].shape[0])
        # self.add_data(Cells=list(df_umi.loc[df_umi['mark'] == 'CB', 'UMI']))
        # self.add_data(UB_num=df_umi[df_umi['mark'] == 'UB'].shape[0])
        # self.add_data(Background=list(df_umi.loc[df_umi['mark'] == 'UB', 'UMI']))
        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))

        title = 'Clonetypes'
        raw_clonotypes = pd.read_csv(self.clonotable, sep=',', index_col=None)
        raw_clonotypes['ClonotypeID'] = raw_clonotypes['clonotype_id'].apply(
            lambda x: x.strip('clonetype'))
        raw_clonotypes['Frequency'] = raw_clonotypes['frequency']
        raw_clonotypes['Proportion'] = raw_clonotypes['proportion'].apply(
            lambda x: f'{round(x*100, 2)}%')
        raw_clonotypes['CDR3_aa'] = raw_clonotypes['cdr3s_aa'].apply(
            lambda x: x.replace(';', '<br>'))
        table_dict = self.get_table_dict(
            title=title,
            table_id='clonetypes',
            df_table=raw_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']]
        )
        self.add_data(table_dict=table_dict)


def summarize(args):

    with Summarize(args, display_title="Cells") as runner:
        runner.run()

def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR',
                        choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument(
            '--barcode_dic', help='10X barcode correspond sgr barcode', required=True)
    return parser
