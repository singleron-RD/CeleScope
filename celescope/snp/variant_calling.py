
from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.tools.step import Step, s_common
from celescope.rna.mkref import Mkref_rna


class Variant_calling(Step):
    """
    ## Features
    - Perform variant calling at single cell level.

    ## Output
    - `{sample}_raw.vcf` Variants are called with bcftools default settings.
    - `{sample}_norm.vcf` Indels are left-aligned and normalized. See https://samtools.github.io/bcftools/bcftools.html#norm for more details.
    """

    def __init__(self, args):
        Step.__init__(self, args)

        # set
        self.barcodes, _num = utils.get_barcode_from_match_dir(args.match_dir)
        self.fasta = Mkref_rna.parse_genomeDir(args.genomeDir)['fasta']
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
        self.debug_subprocess_call(cmd)

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
        """
        max depth 100M
        """
        cmd = (
            f'bcftools mpileup '
            f'-f {self.fasta} '
            f'--threads {self.thread} '
            f'--annotate DP,AD -d 100000000 '
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

    def bcftools_norm(self):
        cmd = (
            'bcftools norm '
            '-m- '
            f'-f {self.fasta} '
            f'{self.raw_vcf_file} '
            '| bcftools norm '
            '-d both '
            f'-o {self.norm_vcf_file} '
        )
        self.debug_subprocess_call(cmd)

    def run(self):

        self.SplitNCigarReads()
        self.call_variants()
        self.bcftools_norm()

    


@utils.add_log
def variant_calling(args):

    with Variant_calling(args) as runner:
        runner.run()


def get_opts_variant_calling(parser, sub_program):

    parser.add_argument("--genomeDir", help=HELP_DICT['genomeDir'], required=True)
    parser.add_argument("--panel", help=HELP_DICT['panel'])
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
