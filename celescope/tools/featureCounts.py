
import os
import re
import pathlib
from tenacity import retry_if_result,retry
from collections import defaultdict


from celescope.rna.mkref import Mkref_rna
from celescope.tools.mkref import Mkref
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT



def check_return_info(return_info):
    """
    check fun return
    """
    if return_info > 1:
        return True
    else:
        return False


class FeatureCounts(Step):
    """
    ## Features
    - Assigning uniquely mapped reads to genomic features with FeatureCounts.
    ## Output
    - `{sample}` Numbers of reads assigned to features (or meta-features).
    - `{sample}_summary` Stat info for the overall summrization results, including number of 
    successfully assigned reads and number of reads that failed to be assigned due to 
    various reasons (these reasons are included in the stat info).
    - `{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam` featureCounts output BAM, 
    sorted by coordinates;BAM file contains tags as following(Software Version>=1.1.8):
        - CB cell barcode
        - UB UMI
        - GN gene name
        - GX gene id
    - `{sample}_name_sorted.bam` featureCounts output BAM, sorted by read name.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # set
        self.gtf = Mkref_rna.parse_genomeDir(self.args.genomeDir)['gtf']
        self.genome = Mkref.parse_genomeDir(self.args.genomeDir)
        self.featureCounts_param = args.featureCounts_param

        #gtf_type
        self.gtf_type = ['exon','gene']

        #stats
        self.metrics_numbers = defaultdict(dict)

    def format_stat(self,gtf_type,featureCount_log_file):
        metrics_strs = ['Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity']
        
        metrics_compiled = {}

        for metrics_str in metrics_strs:
            raw_str = re.escape(metrics_str) + r'.*?(\d+)'
            compiled = re.compile(raw_str, flags=re.S)
            metrics_compiled[metrics_str] = compiled

        with open(featureCount_log_file, 'r') as fh:

            for line in fh:
                line = line.strip()
                if not line:
                    continue

                for metrics_str in metrics_compiled:
                    compiled = metrics_compiled[metrics_str]
                    match = compiled.search(line)
                    if match:
                        self.metrics_numbers[gtf_type][metrics_str] = int(match.group(1))
                        break

    @utils.add_log
    def run_featureCounts(self,outdir,gtf_type):
        cmd = (
            'featureCounts '
            '-s 1 '
            f'-a {self.gtf} '
            f'-o {outdir}/{self.out_prefix.split("/")[-1]} '  # not bam
            '-R BAM '
            f'-T {self.thread} '
            f'-t {gtf_type} '
            f'{self.args.input} '
        )
        if self.featureCounts_param:
            cmd += (" " + self.featureCounts_param)
        self.debug_subprocess_call(cmd)

    @retry(retry=retry_if_result(check_return_info))
    def run(self):
        gtf_type = self.gtf_type[0]
        outdir = f'{self.outdir}/tmp/{gtf_type}'
        pathlib.Path(outdir).mkdir(parents=True,exist_ok=True)
        #out files
        name_sorted_bam = f'{outdir}/{self.out_prefix.split("/")[-1]}_name_sorted.bam'
        input_basename = os.path.basename(self.args.input)
        featureCounts_bam = f'{outdir}/{input_basename}.featureCounts.bam'
        featureCount_log_file = f'{outdir}/{self.out_prefix.split("/")[-1]}.summary'

        self.run_featureCounts(outdir,gtf_type)
        samtools_runner = utils.Samtools(
            in_bam=featureCounts_bam,
            out_bam=featureCounts_bam,
            threads=self.thread,
            debug=self.debug
        )
        samtools_runner.add_tag(self.gtf)
        samtools_runner.temp_sam2bam(by='coord')
        samtools_runner.samtools_sort(
            in_file=featureCounts_bam,
            out_file=name_sorted_bam,
            by='name',
        )
        self.format_stat(gtf_type,featureCount_log_file)

        if (judge_cd:=len(self.gtf_type)):
            self.gtf_type.pop(0)
            return judge_cd

    @utils.add_log
    def creat_stat(self):
        """
        creat stat
        """
        metrics_numbers = {}
        metrics_numbers['Assigned_exon'] = self.metrics_numbers['exon']['Assigned']
        metrics_numbers['Assigned_intron'] = self.metrics_numbers['gene']['Assigned']-self.metrics_numbers['exon']['Assigned']
        metrics_numbers['Assigned_intergenic'] = self.metrics_numbers['exon']['Unassigned_NoFeatures']
        metrics_numbers['Unassigned_Ambiguity'] = self.metrics_numbers['exon']['Unassigned_Ambiguity'] - metrics_numbers['Assigned_intron']
        
        total = sum(metrics_numbers.values())
        
        self.add_metric(
            name='Genome',
            value=self.genome['genome_name'],
        )
        self.add_metric(
            name='Feature Type',
            value=self.args.gtf_type.capitalize(),
        )
        self.add_metric(
            name='Reads Assigned To Exonic Regions',
            value=metrics_numbers['Assigned_exon'],
            total=total,
            help_info='Reads that can be successfully assigned to exonic regions'
        )
        self.add_metric(
            name='Reads Assigned To Intronic Regions',
            value=metrics_numbers['Assigned_intron'],
            total=total,
            help_info='Reads that can be successfully assigned to intron regions'
        )
        self.add_metric(
            name='Reads Assigned To Intergenic Regions',
            value=metrics_numbers['Assigned_intergenic'],
            total=total,
            help_info='Reads that can be successfully assigned to intergenic regions'
        )
        self.add_metric(
            name='Reads Unassigned Ambiguity',
            value=metrics_numbers['Unassigned_Ambiguity'],
            total=total,
            help_info='Alignments that overlap two or more features'
        )


    def clean_tmp(self):
        """
        remove tmp dir
        """
        use_gt = self.args.gtf_type
        retain_dir = f'{self.outdir}/tmp/{use_gt}'
        cmd = (f'mv {retain_dir}/* {self.outdir};rm -rf {self.outdir}/tmp')
        self.debug_subprocess_call(cmd)


@utils.add_log
def featureCounts(args):
    with FeatureCounts(args) as runner:
        runner.run()
        runner.creat_stat()
        runner.clean_tmp()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument(
        '--gtf_type',
        help='Specify feature type in GTF annotation',
        default='exon'
    )
    parser.add_argument('--genomeDir', help=HELP_DICT['genomeDir'])
    parser.add_argument('--featureCounts_param', help=HELP_DICT['additional_param'], default="")

    if sub_program:
        parser.add_argument('--input', help='Required. BAM file path.', required=True)
        parser = s_common(parser)
    return parser
