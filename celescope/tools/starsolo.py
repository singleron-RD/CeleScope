import sys
import subprocess

from celescope.tools.__init__ import PATTERN_DICT
from celescope.__init__ import ROOT_PATH, HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.tools.barcode import Chemistry, Barcode
from celescope.tools import utils


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
        
        if chemistry != 'customized':
            self.whitelist_str = " ".join(Chemistry.get_whitelist(chemistry))
            pattern = PATTERN_DICT[chemistry]
        else:
            self.whitelist_str = args.whitelist
            pattern = self.args.pattern
        self.cb_pos, self.umi_pos, self.umi_len = self.get_solo_pos(pattern)

        if args.cell_calling_method == 'EmptyDrops_CR':
            self.cell_filter = 'EmptyDrops_CR'
        elif args.cell_calling_method == 'auto':
            self.cell_filter = 'CellRanger2.2'
        
        # output files
        solo_dir = f'{self.outdir}/{self.sample}_Solo.out/'
        self.raw_matrix = f'{solo_dir}/GeneFull_Ex50pAS/raw'
        self.filtered_matrix = f'{solo_dir}/GeneFull_Ex50pAS/filtered'
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
            '--soloCBmatchWLtype 1MM \\\n'
            '--soloFeatures Gene GeneFull GeneFull_Ex50pAS \\\n'
            '--outSAMattributes NH HI nM AS CR UR CB UB GX GN \\\n'
            '--outSAMtype BAM SortedByCoordinate \\\n'
            f'--soloCellFilter {self.cell_filter} \\\n'
            f'--outFileNamePrefix {self.out_prefix}_ \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--clip3pAdapterSeq {self.args.adapter_3p} \\\n'
            f'--outFilterMatchNmin {self.args.outFilterMatchNmin} \\\n'
        )
        sys.stderr.write(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def gzip_matrix(self):
        cmd = f'gzip {self.raw_matrix}/*; gzip {self.filtered_matrix}/*'
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_starsolo()
        self.gzip_matrix()

def starsolo(args):
    with Starsolo(args) as runner:
        runner.run()


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
        default='AAAAAAAAAA',
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
        '--cell_calling_method',
        help=HELP_DICT['cell_calling_method'],
        choices=['auto', 'EmptyDrops_CR'],
        default='EmptyDrops_CR',
    )
    parser.add_argument(
        '--starMem',
        help='Default `30`. Maximum memory that STAR can use.',
        default=30
    )
    if sub_program:
        parser.add_argument('--fq1', help='R1 fastq file. Multiple files are separated by comma.', required=True)
        parser.add_argument('--fq2', help='R2 fastq file. Multiple files are separated by comma.', required=True)
        parser = s_common(parser)

    return parser