import subprocess

import celescope.tools.utils as utils
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.rna.mkref import parse_genomeDir_rna


class Variant_calling(Step):
    """
    Features
    - Perform variant calling at single cell level.

    Output

    `{sample}_VID.tsv` A unique numeric ID is assigned for each variant, 

    - `RID`: Target region ID. This column will be added when `--panel` option were provided.

    `{sample}_CID.tsv` A unique numeric ID is assigned for each cell.

    `{sample}_RID.tsv` A unique numeric ID is assigned for each target region. This file will be created when `--panel` option were provided.


    `{sample}_merged.vcf ` VCF file containing all variants of all cells. `VID` and `CID` are added to the `INFO` column.

    `{sample}_filter.vcf` VCF file after filtering. Invalid `CID`s are removed from the `INFO` column.

    `{sample}_variant_count.tsv`  Reference and variant supporting reads/UMIs counts.

    `{sample}_filter_variant_count.tsv`  Reference and variant supporting reads/UMIs counts after filtering.

    `{sample}_support.mtx` Support matrix in [Matrix Market Exchange Formats](https://math.nist.gov/MatrixMarket/formats.html). Rows 
    are variants(VID) and columns are cells(CID). The value can be 1, 2 or 3.
    
    1 : all reads/UMIs at the position support the ref allele.  
    2 : all reads/UMIs at the position support the alt allele.  
    3 : one or more reads/UMIs support both the alt and the ref allele.  
    """

    def __init__(self, args,display_title=None):
        Step.__init__(self, args,display_title=display_title)

        # set
        self.barcodes, _num = utils.read_barcode_file(args.match_dir)
        self.fasta = parse_genomeDir_rna(args.genomeDir)['fasta']
        self.df_vcf = None
        self.panel = args.panel
        self.bed = utils.get_bed_file_path(self.panel)

        # out
        self.splitN_bam = f'{self.out_prefix}_splitN.bam'
        self.splitN_bam_name_sorted = f'{self.out_prefix}_splitN_name_sorted.bam'

        self.raw_bcf_file = f'{self.out_prefix}_raw.bcf'
        self.raw_vcf_file = f'{self.out_prefix}_raw.vcf'
        self.fixed_header_vcf = f'{self.out_prefix}_fixed.vcf'
        self.norm_vcf_file = f'{self.out_prefix}_norm.vcf'


    @utils.add_log
    def SplitNCigarReads(self):
        cmd = (
            f'gatk '
            f'SplitNCigarReads '
            f'--do-not-fix-overhangs '
            f'-R {self.fasta} '
            f'-I {self.args.bam} '
            f'-O {self.splitN_bam} '
        )
        Variant_calling.SplitNCigarReads.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def fix_header(self):
        cmd = (
            'picard FixVcfHeader '
            f'I={self.raw_vcf_file} '
            f'O={self.fixed_header_vcf} '
        )
        self.debug_subprocess_call(cmd)

    @utils.add_log
    def gatk_norm(self):
        cmd = (
            f'gatk LeftAlignAndTrimVariants '
            f'-R {self.fasta} '
            f'-V {self.fixed_header_vcf} '
            f'-O {self.norm_vcf_file} ' 
            '--split-multi-allelics '
        )
        self.debug_subprocess_call(cmd)

    @utils.add_log
    def call_variants(self):
        cmd = (
            f'bcftools mpileup '
            f'-f {self.fasta} '
            f'--threads {self.thread} '
            f'--annotate DP,AD -d 1000000 '
            f'-o {self.raw_bcf_file} '
            f'{self.splitN_bam} '
        )
        if self.bed:
            cmd += f' --regions-file {self.bed} '
        self.debug_subprocess_call(cmd)

        cmd = (
            f'bcftools call '
            f'-mv -Ov '
            f'-o {self.raw_vcf_file} '
            f'{self.raw_bcf_file} '
        )
        self.debug_subprocess_call(cmd)



    def run(self):

        self.SplitNCigarReads()
        self.call_variants()
        self.fix_header()
        self.gatk_norm()
        self.clean_up()


@utils.add_log
def variant_calling(args):

    with Variant_calling(args,display_title="Variant Calling") as runner:
        runner.run()


def get_opts_variant_calling(parser, sub_program):

    parser.add_argument("--genomeDir", help=HELP_DICT['genomeDir'], required=True)
    parser.add_argument("--panel", help = HELP_DICT['panel'])
    if sub_program:
        parser.add_argument(
            "--bam",
            help='Input BAM file from step `target_metrics`. ',
            required=True
        )
        parser.add_argument(
            "--match_dir",
            help=HELP_DICT['match_dir'],
            required=True
        )
        s_common(parser)
