import os
import json
from collections import defaultdict

import pandas as pd
import pysam

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.plotly_plot import Bar_plot


class Summarize(Step):
    """
    ##Features

    - Summarize contig and clonetypes infomation.

    - Generate Mapped reads and Cells metrics, barcode rank plot, clonotypes table and clonotypes frequency barplot in html.

    Output
    - `filtered_contig_annotations.csv` High-level annotations of each high-confidence, cellular contig.

    - `filtered_contig.fasta` High-confidence contig sequences in cell barcodes.

    - `clonotypes.csv` High-level descriptions of each clonotype.

    - `match_contigs.csv` Consider barcodes match scRNA-Seq library in filtered_contig_annotations.csv.

    - `match_contig.fasta` Consider barcodes match scRNA-Seq library in filtered_contig.fasta.
    
    - `match_clonotypes.csv` Consider barcodes match scRNA-Seq library in clonotypes.csv.

    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.seqtype = args.seqtype
        self.not_split_R2 = args.not_split_R2
        self.assemble_out = args.assemble_out
        self.annotation_out = args.annotation_out
        self.match_out = args.match_out

        with open(args.barcode_convert_json, 'r') as f:
            self.tenX_sgr = json.load(f)

        self.all_bam = f'{self.assemble_out}/all_contig.bam'
        self.clono_table = f'{self.assemble_out}/clonotypes.csv'
        self.filter_contig = pd.read_csv(f'{self.annotation_out}/filtered_contig_annotations.csv')
        self.match_contig = pd.read_csv(f'{self.match_out}/match_contigs.csv')
        self.filter_contig_fasta = f'{self.annotation_out}/filtered_contig.fasta'
        self.match_contig_fasta = f'{self.match_out}/match_contig.fasta'

        if self.seqtype == 'TCR':
            self.chains = ['TRA', 'TRB']
            self.pair = ['TRA_TRB']
        elif self.seqtype == 'BCR':
            self.chains = ['IGH', 'IGL', 'IGK']
            self.pair = ['IGH_IGL', 'IGH_IGK']

        # out
        self.productive = f'{self.outdir}/productive'
    
    @staticmethod
    def get_productive_fasta(in_fasta, out_fasta, productive_contig_id):
        """Generate fasta file that only includes productive contig.

        :param in_fasta: filtered contig fasta file.
        :param out_fasta: output fasta file of productive contig.
        :param productive_contig_id: productive contig id.
        """
        handle_file = open(out_fasta, 'w')
        with pysam.FastxFile(in_fasta, 'r') as f:
            for read in f:
                if read.name in productive_contig_id:
                    handle_file.write(">" + read.name + "\n" + read.sequence + "\n")
        handle_file.close()  

    @staticmethod
    def get_productive_result(filter_contig, match_contig, filter_contig_fasta, match_contig_fasta, outdir):
        """Generate annotation file that only includes productive contig.

        :param filter_contig: filter contig annotation.
        :param match_contig: filter match contig annotation.
        :param filter_contig_fasta: filter contig fasta.
        :param match_contig_fasta: filter match contig fasta
        :param outdir: output path.
        """
        productive_contig = filter_contig[filter_contig['productive'] == True]
        productive_contig.to_csv(f'{outdir}/productive_contig_annotations.csv', sep=',', index=False)

        productive_match_contig = match_contig[match_contig['productive'] == True]
        productive_match_contig.to_csv(f'{outdir}/productive_match_contig_annotations.csv', sep=',', index=False)

        productive_contig_id = set(productive_contig['contig_id'])
        productive_fasta, producitve_match_fasta = f'{outdir}/productive_contig.fasta', f'{outdir}/productive_match_contig.fasta'
        Summarize.get_productive_fasta(filter_contig_fasta, productive_fasta, productive_contig_id)
        Summarize.get_productive_fasta(match_contig_fasta, producitve_match_fasta, productive_contig_id)

    @utils.add_log
    def gen_productive_res(self):
        """Generate productive chain result"""

        os.system(f"mkdir -p {self.productive}")
        self.get_productive_result(self.filter_contig, self.match_contig, self.filter_contig_fasta, self.match_contig_fasta, self.productive)

    @utils.add_log
    def gen_report(self):
        """Summarize in html"""

        stat_dict = pd.read_csv(f'{self.outdir}/../01.barcode/stat.txt', sep=':', header=None)
        read_count = int(stat_dict.iloc[0, 1].replace(',', ''))
        sum_dict = pd.read_csv(f'{self.assemble_out}/metrics_summary.csv', sep=',', index_col=None)
        sum_dict = sum_dict.T.to_dict()
        total_reads = int(sum_dict[0]["Number of Read Pairs"].replace(',', ''))
        cell_nums = len(set(self.filter_contig.barcode))

        _index = 200
        if self.not_split_R2:
            _index = 100

        self.add_metric(
            name='Estimated Number of Cells',
            value=cell_nums,
            help_info="Cells with at least one productive TRA or TRB chain"
        )

        self.add_metric(
            name='Reads Mapped to Any V(D)J Gene',
            value=int(total_reads * (float(sum_dict[0]['Reads Mapped to Any V(D)J Gene'].strip('%'))/_index)),
            total=read_count,
            help_info="Reads mapped to any TCR or BCR genes."
        )

        for chain in self.chains:
            self.add_metric(
                name=f'Reads Mapped to {chain}',
                value=int(total_reads * (float(sum_dict[0][f'Reads Mapped to {chain}'].strip('%'))/_index)),
                total=read_count,
                help_info=f"Reads mapped to {chain} chain. "
            )

        self.add_metric(
            name='Fraction Reads in Cells',
            value=int(total_reads * (float(sum_dict[0]['Fraction Reads in Cells'].strip('%'))/_index)),
            total=read_count,
            help_info="Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes"
        )

        for chain in self.chains:
            mid = self.filter_contig[self.filter_contig['chain']== chain]['umis'].median()
            if mid == mid:
                self.add_metric(
                    name=f'Median used {chain} UMIs per Cell',
                    value=int(mid),
                    help_info=f"Median number of UMIs assigned to a {chain} contig per cell."
                )
            else:
                self.add_metric(
                    name=f'Median used {chain} UMIs per Cell',
                    value=0,
                    help_info=f"Median number of UMIs assigned to a {chain} contig per cell."
                )
    
    @utils.add_log
    def get_plot(self):
        """Barcode rank plot and Clonotypes bar plot"""
        
        all_bam = pysam.AlignmentFile(self.all_bam)
        dic_umi = defaultdict(set)

        for read in all_bam:
            cb = read.get_tag('CB')
            umi = read.get_tag('UB')
            new_cb = utils.reverse_complement(self.tenX_sgr[cb.split('-')[0]])
            dic_umi[new_cb].add(umi)

        df_umi = pd.DataFrame()
        df_umi['barcode'] = list(dic_umi.keys())
        df_umi['UMI'] = [len(dic_umi[i]) for i in dic_umi]
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        sgr_cbs = set(self.filter_contig['barcode'].tolist())
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in sgr_cbs else 'UB')

        df_umi.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(f'{self.outdir}/count.txt'))

        title = 'Clonetypes'
        raw_clonotypes = pd.read_csv(self.clono_table, sep=',', index_col=None)
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
        os.system(f'cp {self.annotation_out}/* {self.outdir}')
        os.system(f'cp {self.match_out}/* {self.outdir}')
        self.gen_productive_res()
        self.gen_report()
        self.get_plot()


def summarize(args):
    with Summarize(args, display_title="Cells") as runner:
        runner.run()


def get_opts_summarize(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--not_split_R2', help='not split R2 reads')
    if sub_program:
        s_common(parser)
        parser.add_argument('--barcode_convert_json', help='json file', required=True)
        parser.add_argument('--assemble_out', help='assemble result', required=True)
        parser.add_argument('--annotation_out', help='annotation result', required=True)
        parser.add_argument('--match_out', help='match result', required=True)
    return parser
