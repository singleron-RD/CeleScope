
import os
import re
import subprocess

from celescope.rna.mkref import parse_genomeDir_rna
from celescope.tools.step import Step, s_common
import celescope.tools.utils as utils


class FeatureCounts(Step):
    """
    Features

    - Assigning uniquely mapped reads to genomic features with FeatureCounts.

    Output
    - `{sample}` Numbers of reads assigned to features (or meta-features).

    - `{sample}_summary` Stat info for the overall summrization results, including number of 
    successfully assigned reads and number of reads that failed to be assigned due to 
    various reasons (these reasons are included in the stat info).

    - `{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam` featureCounts output BAM, 
    sorted by coordinates；BAM file contains tags as following(Software Version>=1.1.8):
        - CB cell barcode
        - UB UMI
        - GN gene name
        - GX gene id

    - `{sample}_name_sorted.bam` featureCounts output BAM, sorted by read name.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # set
        self.gtf = parse_genomeDir_rna(self.args.genomeDir)['gtf']
        self.featureCounts_param = args.featureCounts_param

        # out files
        input_basename = os.path.basename(self.args.input)
        self.featureCounts_bam = f'{self.outdir}/{input_basename}.featureCounts.bam'
        self.name_sorted_bam = f'{self.out_prefix}_name_sorted.bam'
        self.featureCount_log_file = f'{self.out_prefix}.summary'

    def format_stat(self):
        tmp_arr = []
        fh = open(self.featureCount_log_file, 'r')
        with open(self.stat_file, 'w') as stat_fh:
            p1 = re.compile(r'Assigned.*?(\d+)', flags=re.S)
            p2 = re.compile(r'Unassigned_NoFeatures.*?(\d+)', flags=re.S)
            p3 = re.compile(r'Unassigned_Ambiguity.*?(\d+)', flags=re.S)
            for line in fh:
                if line.strip() == '':
                    continue

                m1 = p1.search(line.strip())
                if m1:
                    tmp_arr.append(int(m1.group(1)))

                m2 = p2.search(line)
                if m2:
                    tmp_arr.append(int(m2.group(1)))

                m3 = p3.search(line)
                if m3:
                    tmp_arr.append(int(m3.group(1)))

            total = sum(tmp_arr)
            tmp_arr = [
                '%s(%.2f%%)' %
                (utils.format_number(n), (n + 0.0) / total * 100) for n in tmp_arr]
            for t, s in zip(['Assigned', 'Unassigned_NoFeatures',
                            'Unassigned_Ambiguity'], tmp_arr):
                stat_fh.write('%s: %s\n' % (t, s))
        fh.close()

    @utils.add_log
    def run_featureCounts(self):
        cmd = (
            'featureCounts '
            '-s 1 '
            f'-a {self.gtf} '
            f'-o {self.out_prefix} '  # not bam
            '-R BAM '
            f'-T {self.thread} '
            f'-t {self.args.gtf_type} '
            f'{self.args.input} '
        )
        if self.featureCounts_param:
            cmd += (" " + self.featureCounts_param)
        FeatureCounts.run_featureCounts.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


    def run(self):
        self.run_featureCounts()
        samtools_runner = utils.Samtools(
            in_bam=self.featureCounts_bam,
            out_bam=self.featureCounts_bam,
            threads=self.thread,
            )
        samtools_runner.add_tag(self.gtf)
        samtools_runner.temp_sam2bam(by='name')
        self.format_stat()
        self.clean_up()


@utils.add_log
def featureCounts(args):
    step_name = "featureCounts"
    featureCounts_obj = FeatureCounts(args, step_name)
    featureCounts_obj.run()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument(
        '--gtf_type',
        help='Specify feature type in GTF annotation',
        default='exon'
    )
    parser.add_argument('--genomeDir', help='Required. Genome directory.')
    parser.add_argument('--featureCounts_param', help='Other featureCounts parameters', default="")

    if sub_program:
        parser.add_argument('--input', help='Required. BAM file path.', required=True)
        parser = s_common(parser)
    return parser
