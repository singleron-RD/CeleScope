import os
import subprocess
import pandas as pd 
from celescope.tools import utils
from celescope.tools.step import Step, s_common

__SUB_STEPS__ = ['sequencing', 'mapping', 'cells', 'targeting']

"""
read1, barcode, read2, sample index are associated with R1, R2, R3, I1 respectively.
R1: Read 1
R2: Dual index i5 read(10x Barcode)
R3: Read 2
"""


def get_opts_atac(parser, sub_program):
    parser.add_argument('--ref_path', help='reference path for cellranger-atac', required=True)
    parser.add_argument('--soft_path', help='soft path for cellranger-atac', required=True)
    parser.add_argument('--other_param', help='Other cellranger-atac count parameters.', default="")
    parser.add_argument('--mem', help='memory(G)', default=30)
    if sub_program:
        s_common(parser)
        parser.add_argument('--fqs_dir', help='fastq dir after convert', required=True)
    return parser


def format_value(value):
    """
    >>> value = 0.8941
    >>> format_value(value)
    89.41%
    """
    return str(round(value * 100, 2)) + '%'


class ATAC(Step):

    """
    ##Features
    - Run cellranger-atac
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fqs_dir = os.path.abspath(args.fqs_dir)
        self.mem = args.mem
        self.ref_path = args.ref_path
        self.soft_path = args.soft_path

        # out files
        self.cmd_line = f'{self.outdir}/{self.sample}_cmd_line'

    @utils.add_log
    def cr_atac(self):
        """cellranger-atac count"""
        cmd = (
            f'{self.soft_path} count '
            f'--id={self.sample} '
            f'--reference={self.ref_path} '
            f'--fastqs={self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )
        if self.args.other_param:
            cmd += (self.args.other_param)

        self.cr_atac.logger.info(cmd)
        with open(self.cmd_line, 'w') as f:
            f.write(cmd)
        
        cwd = os.getcwd()
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid can not find '02.atac/stat.txt' error
        os.chdir(cwd)

    def run(self):
        self.cr_atac()


def atac(args):
    with ATAC(args) as runner:
        runner.run()

    with Sequencing(args) as runner:
        runner.run()

    with Mapping(args) as runner:
        runner.run()

    with Cells(args) as runner:
        runner.run()

    with Targeting(args) as runner:
        runner.run()


class cellranger_metrics(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        cellranger_metrics_csv = f'{self.outdir}/{self.sample}/outs/summary.csv'
        df = pd.read_csv(cellranger_metrics_csv, index_col=None)
        self.metrics_dict = df.T.to_dict()[0]

        self.species_list = ['']
        if 'GRCh38-and-mm10' in args.ref_path:
            self.species_list = [' (GRCh38)', ' (mm10)']
    
    def run(self):
        self.add_metrics()


class Sequencing(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
    
    def add_metrics(self):

        name = 'Sequenced read pairs'
        self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info="Total number of sequenced read pairs assigned to the sample."
        )

        name = 'Valid barcodes'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of read pairs with barcodes that match the whitelist after error correction."
        )
        
        name = 'Q30 bases in barcode'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of barcode read (i2) bases with Q-score >= 30."
        )

        name = 'Q30 bases in read 1'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of read 1 bases with Q-score >= 30."
        )

        name = 'Q30 bases in read 2'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of read 2 bases with Q-score >= 30."
        )


class Mapping(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        
    def add_metrics(self):

        name = 'Confidently mapped read pairs'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of sequenced read pairs with mapping quality > 30."
        )

        name = 'Unmapped read pairs'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of sequenced read pairs that have a valid barcode but could not be mapped to the genome."
        )

        name = 'Non-nuclear read pairs'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of sequenced read pairs that have a valid barcode and map to non-nuclear genome contigs, including mitochondria,with mapping quality > 30."
        )

        name = 'Fragments in nucleosome-free regions'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of high-quality fragments smaller than 124 basepairs."
        )

        name = 'Fragments flanking a single nucleosome'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of high-quality fragments between 124 and 296 basepairs."
        )


class Cells(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
    
    def add_metrics(self):
        for species in self.species_list:
            name = f'Estimated number of cells{species}'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info="The total number of barcodes identified as cells."
            )

            name = f'Mean raw read pairs per cell{species}'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info="Total number of read pairs divided by the number of cell barcodes"
            )

            name = f'Median high-quality fragments per cell{species}'
            self.add_metric(
                name=name,
                value=self.metrics_dict[name],
                help_info="The median number of high-quality fragments per cell barcode."
            )

        name = 'Fraction of high-quality fragments in cells'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of high-quality fragments with a valid barcode that are associated with cell-containing partitions. High-quality fragments are defined as read pairs with a valid barcode that map to the nuclear genome with mapping quality > 30, are not chimeric and not duplicate."
        )

        name = 'Fraction of transposition events in peaks in cells'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of transposition events that are associated with cell-containing partitions and fall within peaks. Transposition events are located at both ends of all high-quality fragments. This metric measures the percentage of such events that overlap with peaks."
        )


class Targeting(cellranger_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    def add_metrics(self):
        name = 'Number of peaks'
        self.add_metric(
            name=name,
            value=self.metrics_dict[name],
            help_info="Total number of peaks on primary contigs either detected by the pipeline or input by the user."
        )
        
        name = 'Fraction of genome in peaks'
        self.add_metric(
            name=name,
            value=format_value(self.metrics_dict[name]),
            help_info="Fraction of bases in primary contigs that are defined as peaks."
        )
        
        name = 'TSS enrichment score'
        self.add_metric(
            name=name,
            value=round(self.metrics_dict[name], 2),
            help_info="Maximum value of the transcription-start-site (TSS) profile.The TSS profile is the summed accessibility signal (defined as number of cut sites per base) in a window of 2,000 bases around all the annotated TSSs, normalized by the minimum signal in the window."
        )

        for species in self.species_list:
            name = f'Fraction of high-quality fragments overlapping TSS{species}'
            self.add_metric(
                name=name,
                value=format_value(self.metrics_dict[name]),
                help_info="Fraction of high-quality fragments in cell barcodes that overlap transcription start sites (TSS)."
            )

            name = f'Fraction of high-quality fragments overlapping peaks{species}'
            self.add_metric(
                name=name,
                value=format_value(self.metrics_dict[name]),
                help_info="Fraction of high-quality fragments in cell barcodes that overlap called peaks."
            )