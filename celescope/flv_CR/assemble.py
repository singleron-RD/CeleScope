import os
import subprocess
import pysam
import pandas as pd
import json

from collections import defaultdict

from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.flv_trust4.__init__ import CHAIN, PAIRED_CHAIN


__SUB_STEPS__ = ['mapping', 'cells', 'annotation']


class Assemble(Step):
    """
    ## Features

    - TCR/BCR Assemble by Cellranger.

    - Generate Mapping, Cells, V(D)J annotations metrics in html.

    ## Output
    
    - `03.assemble/{sample}` Cellranger vdj results.

    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fqs_dir = os.path.abspath(args.fqs_dir)
        self.mem = args.mem
        self.other_param = args.other_param
        self.ref_path = args.ref_path
        self.soft_path = args.soft_path

        # out files
        self.cmd_line = f'{self.outdir}/{self.sample}_cmd_line'


    @utils.add_log
    def assemble(self):
        """Cellranger vdj"""
        cmd = (
            f'{self.soft_path} vdj '
            f'--id={self.sample} '
            f'--reference={self.ref_path} '
            f'--fastqs={self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )

        if self.other_param:
            cmd += (" " + self.other_param)

        self.assemble.logger.info(cmd)
        with open(self.cmd_line, 'w') as f:
            f.write(cmd)
        
        cwd = os.getcwd()
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid can not find '03.assemble/stat.txt' error
        os.chdir(cwd)

    def run(self):
        self.assemble()

def assemble(args):
    with Assemble(args) as runner:
        runner.run()

    with Mapping(args) as runner:
        runner.run()

    with Cells(args) as runner:
        runner.run()

    with Annotation(args) as runner:
        runner.run()


class cellranger_metrics(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        cellranger_metrics_csv = f'{self.outdir}/{self.sample}/outs/metrics_summary.csv'
        df = pd.read_csv(cellranger_metrics_csv, index_col=None)
        self.metrics_dict = df.T.to_dict()[0]

        self.seqtype = args.seqtype
        self.chains = CHAIN[self.seqtype]
        self.pairs = PAIRED_CHAIN[self.seqtype]

    def run(self):
        self.add_metrics()


class Mapping(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)


    def add_metrics(self):
        name = 'Reads Mapped to Any V(D)J Gene'
        self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info="Fraction of reads that partially or wholly map to any germline V(D)J gene segment."
        )

        for chain in self.chains:
            name = f'Reads Mapped to {chain}'
            self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info=f"Fraction of reads that map partially or wholly to a germline {chain} gene segment."
        )


class Cells(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.barcode_convert_json = f'{args.fqs_dir}/barcode_convert.json'

        with open(self.barcode_convert_json, 'r') as f:
            self.tenX_sgr = json.load(f)

        self.all_bam = f'{self.outdir}/{self.sample}/outs/all_contig.bam'
        self.filter_contig_file = pd.read_csv(f'{self.outdir}/{self.sample}/outs/filtered_contig_annotations.csv')
        
        # out
        self.count_file = f'{self.outdir}/count.txt'

    def add_metrics(self):
        name = 'Estimated Number of Cells'
        self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info="The number of barcodes estimated to be associated with cells that express targeted V(D)J transcripts."
        )

        name = 'Fraction Reads in Cells'
        self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info="Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes."
        )

        name = 'Mean Reads per Cell'
        try:
            value = self.metrics_dict[name]
        except KeyError:
            value = self.metrics_dict['Mean Read Pairs per Cell']
        self.add_metric(
            name=name,
            value=value,
            help_info="Number of input read pairs divided by the estimated number of cells."
        )

        name = 'Mean Used Reads per Cell'
        try:
            value = self.metrics_dict[name]
        except KeyError:
            value = self.metrics_dict['Mean Used Read Pairs per Cell']
        self.add_metric(
            name=name,
            value=value,
            help_info="Mean number of read pairs used in assembly per cell-associated barcode."
        )

        for chain in self.chains:
            median_value = self.filter_contig_file[self.filter_contig_file['chain']==chain]['umis'].median()
            if median_value == median_value:
                self.add_metric(
                    name=f'Median used {chain} UMIs per Cell',
                    value=int(median_value),
                    help_info=f"Median number of UMIs assigned to a {chain} contig per cell."
                )
            else:
                self.add_metric(
                    name=f'Median used {chain} UMIs per Cell',
                    value=0,
                    help_info=f"Median number of UMIs assigned to a {chain} contig per cell."
                )

    @utils.add_log
    def Barcode_rank_plot(self):
        dic_umi = defaultdict(set)
        
        with pysam.AlignmentFile(self.all_bam) as fh:
            for read in fh:
                cb = read.get_tag('CB')
                umi = read.get_tag('UB')
                dic_umi[cb].add(umi)

        df_umi = pd.DataFrame()
        df_umi['barcode'] = list(dic_umi.keys())
        df_umi['UMI'] = [len(dic_umi[i]) for i in dic_umi]
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        cbs = set(self.filter_contig_file['barcode'])
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in cbs else 'UB')
        df_umi['barcode'] = df_umi['barcode'].apply(lambda x : self.tenX_sgr[x.split('-')[0]])

        df_umi.to_csv(self.count_file, sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.count_file))

    def run(self):
        self.add_metrics()
        self.Barcode_rank_plot()


class Annotation(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
    
    def add_metrics(self):
        name = 'Cells With Productive V-J Spanning Pair'
        self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info="Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair.A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region."
        )

        for pair in self.pairs:
            chain1, chain2 = pair.split('_')[0], pair.split('_')[1]
            name = f'Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info = f"Fraction of cell-associated barcodes with at least one productive contig for each chain of the ({chain1}, {chain2}) receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region."
            )
        
        for chain in self.chains:
            name = f'Cells With {chain} Contig'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info=f"Fraction of cell-associated barcodes with at least one {chain} contig annotated as a full or partial V(D)J gene."
            )

            name = f'Cells With CDR3-annotated {chain} Contig'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info=f"Fraction of cell-associated barcodes with at least one {chain} contig where a CDR3 was detected."
            )

            name = f'Cells With V-J Spanning {chain} Contig'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info=f"Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to the 3' end of the J region for {chain}."
            )

            name = f'Cells With Productive {chain} Contig'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info=f"Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for {chain}, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region."
            )

    
def get_opts_assemble(parser, sub_program):
    parser.add_argument('--ref_path', help='reference path for cellranger')
    parser.add_argument('--soft_path', help='soft path for cellranger')
    parser.add_argument('--other_param', help='Other cellranger parameters.', default="")
    parser.add_argument('--mem', help='memory(G)', default=10)
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir after convert', required=True)
    return parser

