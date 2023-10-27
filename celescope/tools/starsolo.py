import sys
import subprocess
from collections import Counter

import pandas as pd
import pysam
from scipy import interpolate
import numpy as np
import math

from celescope.tools.__init__ import PATTERN_DICT, FILTERED_MATRIX_DIR_SUFFIX, COUNTS_FILE_NAME 
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.tools.barcode import Chemistry, Barcode
from celescope.tools import utils
from celescope.tools.make_ref import MakeRef
from celescope.tools.matrix import CountMatrix
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.cells import Cells_metrics

SAM_attributes = 'NH HI nM AS CR UR CB UB GX GN '
MIN_CELL = 500
MAX_CELL = 60000
MIN_BEAD = 100000
MAX_BEAD = 300000

class Starsolo(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        fq1_list = args.fq1.split(",")
        fq2_list = args.fq2.split(",")
        fq1_number = len(fq1_list)
        fq2_number = len(fq2_list)
        if fq1_number != fq2_number:
            sys.exit('fastq1 and fastq2 do not have same file number!')
        self.read_command = 'cat'
        if str(fq1_list[0]).endswith('.gz'):
            self.read_command = 'zcat'
        if args.chemistry == 'auto':
            ch = Chemistry(args.fq1)
            chemistry_list = ch.check_chemistry()
            if len(set(chemistry_list)) != 1:
                sys.exit('multiple chemistry found!' + str(chemistry_list))
            chemistry = chemistry_list[0]
        else:
            chemistry = args.chemistry
        self.chemistry = chemistry
        
        if chemistry != 'customized':
            self.whitelist_str = " ".join(Chemistry.get_whitelist(chemistry))
            pattern = PATTERN_DICT[chemistry]
        else:
            self.whitelist_str = args.whitelist
            pattern = self.args.pattern
        self.cb_pos, self.umi_pos, self.umi_len = self.get_solo_pos(pattern)
        self.pattern = pattern
       
        # output files
        self.solo_out_dir = f'{self.outdir}/{self.sample}_Solo.out/'
        solo_dir = f'{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS'
        self.raw_matrix = f'{solo_dir}/raw'
        self.filtered_matrix = f'{solo_dir}/filtered'
        self.summary_file = f'{solo_dir}/Summary.csv'
        bam = f'{self.outdir}/{self.sample}_Aligned.sortedByCoord.out.bam'

        # outs
        self.outs = [self.raw_matrix, self.filtered_matrix, bam]
        
    
    @staticmethod
    def get_solo_pos(pattern):
        # returns: cb_pos, umi_pos, umi_len
        pattern_dict = Barcode.parse_pattern(pattern)
        cb_pos = ' '.join([f'0_{l}_0_{r-1}' for l, r in pattern_dict["C"]])
        if len(pattern_dict['U']) != 1:
            sys.exit(f'Error: Wrong pattern:{pattern}. \n Solution: fix pattern so that UMI only have 1 position.\n')
        ul, ur = pattern_dict["U"][0]
        umi_pos = f'0_{ul}_0_{ur-1}'
        umi_len = ur - ul
        return cb_pos, umi_pos, umi_len

    def run_starsolo(self):
        sa = SAM_attributes + self.args.SAM_attributes
        cmd = (
            'STAR \\\n'
            f'--genomeDir {self.args.genomeDir} \\\n'
            f'--readFilesIn {self.args.fq2} {self.args.fq1} \\\n'
            f'--readFilesCommand {self.read_command} \\\n'
            f'--soloCBwhitelist {self.whitelist_str} \\\n'
            f'--soloType CB_UMI_Complex \\\n'
            f'--soloCBposition {self.cb_pos} \\\n'
            f'--soloUMIposition {self.umi_pos} \\\n'
            f'--soloUMIlen {self.umi_len} \\\n'
            f'--soloCellFilter {self.args.soloCellFilter} \\\n'
            f'--outFileNamePrefix {self.out_prefix}_ \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--clip3pAdapterSeq {self.args.adapter_3p} \\\n'
            f'--outFilterMatchNmin {self.args.outFilterMatchNmin} \\\n'
            f'--soloFeatures {self.args.soloFeatures} \\\n'
            f'--outSAMattributes {sa} \\\n'
            '--soloCBmatchWLtype 1MM \\\n'
            '--outSAMtype BAM SortedByCoordinate \\\n'
            '--soloCellReadStats Standard \\\n'
        )
        if self.args.STAR_param:
            cmd += self.args.STAR_param
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)
        cmd = f'chmod -R 755 {self.solo_out_dir}'
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def gzip_matrix(self):
        cmd = f'gzip {self.raw_matrix}/*; gzip {self.filtered_matrix}/*'
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_Q30_cb_UMI(self):
        fq1_list = self.args.fq1.split(",")
        pattern_dict = Barcode.parse_pattern(self.pattern)
        cb_10k, umi_10k, cb_10k_100k, umi_10k_100k = Counter(), Counter(), Counter(), Counter()
        n = 0
        with pysam.FastxFile(fq1_list[0], persist=False) as fq1:
            for entry in fq1:
                n += 1
                if n > 10 ** 5:
                    break
                qual = entry.quality
                cb_qual = Barcode.get_seq_str(qual, pattern_dict['C'])
                umi_qual = Barcode.get_seq_str(qual, pattern_dict['U'])
                if n <= 10 ** 4:
                    cb_10k.update(cb_qual)
                    umi_10k.update(umi_qual)
                else:
                    cb_10k_100k.update(cb_qual)
                    umi_10k_100k.update(umi_qual)
        
        cb_qual_counter = cb_10k
        umi_qual_counter = umi_10k
        if cb_10k_100k:
            cb_qual_counter = cb_10k_100k
            umi_qual_counter = umi_10k_100k
        q30_cb = sum([cb_qual_counter[k] for k in cb_qual_counter if k >= Barcode.ord2chr(
            30)]) / float(sum(cb_qual_counter.values()))
        q30_umi = sum([umi_qual_counter[k] for k in umi_qual_counter if k >= Barcode.ord2chr(
            30)]) / float(sum(umi_qual_counter.values()))
        return q30_cb, q30_umi


    def run(self):
        self.run_starsolo()
        self.gzip_matrix()
        q30_cb, q30_umi = self.get_Q30_cb_UMI()
        return q30_cb, q30_umi, self.chemistry


def starsolo(args):

    with Starsolo(args) as runner:
        q30_cb, q30_umi, chemistry = runner.run()

    with Mapping(args) as runner:
        valid_reads, corrected = runner.run()
    
    with Cells(args) as runner:
        n_reads, q30_RNA = runner.run(chemistry, valid_reads)
    
    with Demultiplexing(args) as runner:
        runner.run(valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA)


class Mapping(Step):
    # only add metrics
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        solo_dir = f'{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS'
        self.cellReadsStats = f'{solo_dir}/CellReads.stats'
        self.filtered_matrix = f'{self.outs_dir}/{FILTERED_MATRIX_DIR_SUFFIX}'
        self.counts_file = f'{solo_dir}/{COUNTS_FILE_NAME }'
        self.genome = MakeRef.get_config(args.genomeDir)['meta']['genome_name']

        self.outs = [self.counts_file]


    @utils.add_log
    def run(self):
        df = pd.read_csv(self.cellReadsStats, sep='\t', header=0, index_col=0)
        df = df.iloc[1:,] # skip first line cb not pass whitelist
        df_count = df.loc[:,['nUMIunique','countedU']] # keep dataframe format
        df_count.rename(columns={'nUMIunique': 'UMI'}, inplace=True)
        df = df.loc[:,['cbMatch', 'cbPerfect','genomeU', 'genomeM', 'exonic', 'intronic','exonicAS','intronicAS','countedU']]
        s= df.sum()
        # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
        valid = int(s['cbMatch'])
        perfect = int(s['cbPerfect'])
        corrected = valid - perfect
        genomeU= int(s['genomeU'])
        genomeM = int(s['genomeM'])
        exonic = int(s['exonic'])
        intronic = int(s['intronic'])
        antisense = int(s['exonicAS'] + s['intronicAS'])
        intergenic = genomeM + genomeU - exonic - intronic - antisense
        countedU = int(s['countedU'])
        del df

        self.add_metric(
            name='genome',
            value=self.genome,
        )

        self.add_metric(
            name='Reads mapped to unique loci',
            value=genomeU / valid,
            value_type='fraction',
            help_info='Reads that mapped uniquely to the genome.'
        )

        self.add_metric(
            name='Reads mapped to multiple loci',
            value=genomeM / valid,
            value_type='fraction',
            help_info='Reads that mapped to multiple loci in the genome'
        )
        unique_transcriptome = countedU / valid
        self.add_metric(
            name='Reads mapped uniquely to Transcriptome',
            value=unique_transcriptome,
            value_type='fraction',
            help_info='Reads that mapped to a unique gene in the transcriptome. These reads are used for UMI counting.'
        )
        self.add_metric(
            name='Reads assigned to exonic regions',
            value=exonic / valid,
            value_type='fraction',
            help_info='Reads that assigned to exonic regions of genes',
        )
        self.add_metric(
            name='Reads assigned to intronic regions',
            value=intronic / valid,
            value_type='fraction',
            help_info='Reads that assigned to intronic regions of genes',
        )
        self.add_metric(
            name='Reads assigned to intergenic regions',
            value=intergenic / valid,
            value_type='fraction',
            help_info='Reads that can not be assigned to a gene will be considered as intergenic reads.',
        )
        self.add_metric(
            name='Reads assigned Antisense to gene',
            value=antisense / valid,
            value_type='fraction',
            help_info='Reads that assigned to the opposite strand of genes',
        )

        df_count.sort_values(by='UMI', ascending=False, inplace=True)
        cbs = CountMatrix.read_barcodes(self.filtered_matrix)
        df_count['mark'] = 'UB'
        for cb in cbs:
            df_count.loc[cb, 'mark'] = 'CB'
        df_count.to_csv(self.counts_file, sep='\t', index=True)

        return valid, corrected


class Cells(Cells_metrics):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)
        solo_dir = f'{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS'
        self.summary_file = f'{solo_dir}/Summary.csv'
        self.counts_file = f'{self.outs_dir}/{COUNTS_FILE_NAME}'


    @utils.add_log
    def parse_summary_add_metrics(self, valid_reads):
        df = pd.read_csv(self.summary_file, index_col=0, header=None)
        s = df.iloc[:,0]
        n_cells = int(s["Estimated Number of Cells"])
        fraction_reads_in_cells = float(s["Fraction of Unique Reads in Cells"])
        mean_used_reads_per_cell = int(s["Mean Reads per Cell"])
        median_umi_per_cell  = int(s["Median UMI per Cell"])
        median_genes_per_cell = int(s["Median GeneFull_Ex50pAS per Cell"])
        total_genes = int(s["Total GeneFull_Ex50pAS Detected"])
        saturation = float(s["Sequencing Saturation"])
        n_reads = int(s["Number of Reads"])
        q30_RNA = float(s["Q30 Bases in RNA read"])

        self.add_cells_metrics(n_cells, fraction_reads_in_cells, mean_used_reads_per_cell, median_umi_per_cell, total_genes, median_genes_per_cell, saturation, valid_reads)
        return n_reads, q30_RNA, fraction_reads_in_cells

    def add_curve_metrics(self, chemistry):
        """
        Returns:
            low barcodes
            loss of cliff
            loss of knee
        """
        df = pd.read_csv(self.counts_file, header=0,sep='\t',index_col=0)
        total = df.shape[0]
        _x = list(range(1, total+1))
        _x = np.log10(_x)
        _y = df['UMI']
        _y = np.log10(_y+1)
        b = interpolate.splrep(_x, _y, s=10)
        d1 = interpolate.splev(_x,b,der=1)
        d2 = interpolate.splev(_x,b,der=2)
        def get_d1_min(x1,x2):
            """
                x1,x2: start and end bc index of the curve
                returns: bc_index, minimum d1
            """
            i = np.argmin(d1[x1:x2+1]) + x1
            return i, d1[i]

        def get_knee_point(x1, x2):
            """
                cliff point(first knee)
                returns: bc_index
            """
            curvature = d2/(1 + d1 ** 2) ** 1.5
            i = np.argmin(curvature[x1:x2+1]) + x1
            return i

        self.add_metric(
            name='Total barcodes',
            value=total,
            help_info='Total deteced barcodes in whitelist.',
            show=False,
        )
        bc_total = math.inf
        if chemistry.startswith('scopeV2'):
            bc_total = 96 ** 3
        elif chemistry.startswith('scopeV3'):
            bc_total = 96 ** 3 * 2
        min_total = bc_total / 2
        low_barcodes = total < min_total
        self.add_metric(
            name='Low total barcodes',
            value=str(low_barcodes),
            show=False,
            help_info='Total number of detected barcodes is lower than expected. This can can be caused by a sample clog.'
        )

        if total < MAX_BEAD:
            sys.stderr.write(f'Warning: Total number of detected barcodes is less than {MAX_BEAD}.\n')
            return True, True, True
        d1_cliff_cell, d1_cliff = get_d1_min(MIN_CELL, MAX_CELL)
        d1_knee_cell, d1_knee = get_d1_min(MIN_BEAD, MAX_BEAD)
        cliff_cell = get_knee_point(MIN_CELL, MAX_CELL)
        self.add_metric(
            name='Minimum derivative of cliff',
            value=round(d1_cliff, 4),
            show=False,
        )
        self.add_metric(
            name='Barcode number at minimum derivative of cliff',
            value=int(d1_cliff_cell),
            show=False,
        )
        self.add_metric(
            name='Barcode number at cliff',
            value=int(cliff_cell),
            show=False,
            help_info='Using the method from emptyDrops to determine the first knee. It is not accurate when the sample shows loss of cliff.'
        )
        loss_of_cliff = d1_cliff > -2
        self.add_metric(
            name='Loss of cliff',
            value=str(loss_of_cliff),
            show=False,
            help_info='If the minimum first derivative around cliff is greater than -2, it is considered as loss of cliff (loss of single cell behavior).'
        )
        self.add_metric(
            name='Minimum derivative of knee',
            value=round(d1_knee, 4),
            show=False,
        )
        self.add_metric(
            name='Barcode number at minimum derivative of knee',
            value=int(d1_knee_cell),
            show=False,
        )
        loss_of_knee = d1_knee > -2
        self.add_metric(
            name='Loss of knee',
            value=str(loss_of_knee),
            show=False,
            help_info='If the minimum first derivative around knee is greater than -2, it is considered as loss of knee.'
        )

        return low_barcodes, loss_of_cliff, loss_of_knee

    @staticmethod
    def get_curve_quality(fraction_reads_in_cells, low_barcodes, loss_of_cliff, loss_of_knee):
        """

        """
        if any([loss_of_cliff, low_barcodes, loss_of_knee]):
            return "fail"
        if fraction_reads_in_cells >= 0.88:
            return 'super'
        elif fraction_reads_in_cells >= 0.75:
            return 'good'
        return 'pass'


    def run(self, chemistry, valid_reads):
        n_reads, q30_RNA, fraction_reads_in_cells = self.parse_summary_add_metrics(valid_reads)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.counts_file))
        low_barcodes, loss_of_cliff, loss_of_knee = self.add_curve_metrics(chemistry)
        quality = self.get_curve_quality(fraction_reads_in_cells, low_barcodes, loss_of_cliff, loss_of_knee)
        self.add_metric(
            name='Barcode rank curve quality', 
            value=quality,
            show=False,
        )
        return n_reads, q30_RNA


class Demultiplexing(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

    def run(self, valid_reads, n_reads, corrected, q30_cb, q30_umi, q30_RNA):
        self.add_metric(
            name='Raw Reads',
            value=n_reads,
            help_info='total reads from FASTQ files'
        )
        self.add_metric(
            name='Valid Reads',
            value=valid_reads / n_reads,
            value_type='fraction',
            help_info='fraction of reads with valid barcode and UMI'
        )

        self.add_metric(
            name='Corrected Barcodes',
            value=corrected / valid_reads,
            value_type='fraction',
            help_info='fraction of valid reads with corrected barcodes. Barcodes are corrected to the whitelist sequence that is 1 Hamming-distance away.'
        )

        self.add_metric(
            name='Q30 of Barcodes',
            value=q30_cb,
            value_type='fraction',
            help_info='Fraction of barcode bases with quality score >= 30',
        )

        self.add_metric(
            name='Q30 of UMI',
            value=q30_umi,
            value_type='fraction',
            help_info='Fraction of UMI bases with quality score >= 30',
        )

        self.add_metric(
            name='Q30 of RNA Reads',
            value=q30_RNA,
            value_type='fraction',
            help_info='Fraction of RNA read bases with quality score >= 30',
        )


def get_opts_starsolo(parser, sub_program=True):
    parser.add_argument(
        '--chemistry',
        help='Predefined (pattern, barcode whitelist, linker whitelist) combinations. ' + HELP_DICT['chemistry'],
        choices=list(PATTERN_DICT.keys()),
        default='auto'
    )
    parser.add_argument(
        '--pattern',
        help="""The pattern of R1 reads, e.g. `C8L16C8L16C8L1U12T18`. The number after the letter represents the number 
        of bases.  
        - `C`: cell barcode  
        - `L`: linker(common sequences)  
        - `U`: UMI    
        - `T`: poly T""",
    )
    parser.add_argument(
        '--whitelist',
        help='Cell barcode whitelist file path, one cell barcode per line.'
    )
    parser.add_argument(
        '--adapter_3p',
        help='Adapter sequence to clip from 3 prime. Multiple sequences are seperated by space',
        default='AAAAAAAAAAAA',
    )
    parser.add_argument(
        '--genomeDir',
        help=HELP_DICT['genomeDir'],
    )
    parser.add_argument(
        '--outFilterMatchNmin',
        help="""Alignment will be output only if the number of matched bases 
is higher than or equal to this value.""",
        default=50,
    )
    parser.add_argument(
        '--soloCellFilter',
        help='The same as the argument in STARsolo',
        default='EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000',
    )
    parser.add_argument(
        '--starMem',
        help='Maximum memory that STAR can use.',
        default=32
    )
    parser.add_argument('--STAR_param', help=HELP_DICT['additional_param'], default="")
    parser.add_argument(
        '--SAM_attributes', 
        help=f'Additional attributes(other than {SAM_attributes}) to be added to SAM file',
        default="")
    parser.add_argument(
        '--soloFeatures', 
        help='The same as the argument in STARsolo',
        default='Gene GeneFull_Ex50pAS',
    )
    if sub_program:
        parser.add_argument('--fq1', help='R1 fastq file. Multiple files are separated by comma.', required=True)
        parser.add_argument('--fq2', help='R2 fastq file. Multiple files are separated by comma.', required=True)
        parser = s_common(parser)

    return parser