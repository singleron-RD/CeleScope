
import os
import pathlib

from collections import defaultdict

import pandas as pd

from celescope.rna.mkref import Mkref_rna
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT



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
        self.featureCounts_param = args.featureCounts_param

        #gtf_type
        self.gtf_types = ['exon','gene']

        #stats
        self.feature_log_dict = defaultdict(dict)


    @staticmethod
    def read_log(log_file):
        """
        Args:
            log_file: featureCounts log summary file
        Returns:
            log_dict: {'Assigned': 123, ...}
        """
        # skip first line
        df = pd.read_csv(log_file, sep='\t', header=None, names=['name', 'value'], skiprows=1)
        log_dict = df.set_index('name')['value'].to_dict()
        return log_dict


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
            '2>&1 '
        )
        if self.featureCounts_param:
            cmd += (" " + self.featureCounts_param)
        self.debug_subprocess_call(cmd)



    def run(self):
        for gtf_type in self.gtf_types:
            outdir = f'{self.outdir}/tmp/{gtf_type}'
            pathlib.Path(outdir).mkdir(parents=True,exist_ok=True)
            #out files
            name_sorted_bam = f'{outdir}/{self.out_prefix.split("/")[-1]}_name_sorted.bam'
            input_basename = os.path.basename(self.args.input)
            featureCounts_bam = f'{outdir}/{input_basename}.featureCounts.bam'
            log_file = f'{outdir}/{self.out_prefix.split("/")[-1]}.summary'

            self.run_featureCounts(outdir,gtf_type)
            if gtf_type == self.args.gtf_type:
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
            self.feature_log_dict[gtf_type] = FeatureCounts.read_log(log_file)


    @utils.add_log
    def add_metrics(self):
        total = sum(self.feature_log_dict['exon'].values())

        Assigned_exon = self.feature_log_dict['exon']['Assigned']
        Assigned_intergenic = self.feature_log_dict['gene']['Unassigned_NoFeatures']
        """
        https://academic.oup.com/nargab/article/2/3/lqaa073/5910008
        Approximately 15% of genes had exon counts that were greater than genebody counts (by a median value of eight counts). This was due to our conservative approach of excluding reads that overlapped features in multiple genes during the read summarization step by featureCounts using the argument allowMultiOverlap=FALSE. Under this strategy, some reads were counted towards the exon count set but not the genebody count set. This happens when a read overlaps the exon in one gene and the intron of another gene—it is counted towards exon counts but not genebody counts due to its overlap of multiple genebodies but not multiple exons.
        """
        Unassigned_ambiguity = self.feature_log_dict['exon']['Unassigned_Ambiguity']
        Assigned_intron = total - Assigned_exon - Assigned_intergenic - Unassigned_ambiguity           
       
        
        self.add_metric(
            name='Feature Type',
            value=self.args.gtf_type.capitalize(),
            help_info='Specified by `--gtf_type`. For snRNA-seq, you need to add `--gtf_type gene` to include reads mapped to intronic regions. Staring from CeleScope v1.12.0, the default value of gtf_type is changed from `exon` to `gene`.'
        )
        self.add_metric(
            name='Reads Assigned To Exonic Regions',
            value=Assigned_exon,
            total=total,
            help_info='Reads that can be successfully assigned to exonic regions'
        )
        self.add_metric(
            name='Reads Assigned To Intronic Regions',
            value=Assigned_intron,
            total=total,
            help_info='Reads that can be successfully assigned to intronic regions'
        )
        self.add_metric(
            name='Reads Assigned To Intergenic Regions',
            value=Assigned_intergenic,
            total=total,
            help_info='Reads that can be successfully assigned to intergenic regions'
        )
        self.add_metric(
            name='Reads Unassigned Ambiguity',
            value=Unassigned_ambiguity,
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
        runner.add_metrics()
        runner.clean_tmp()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument(
        '--gtf_type',
        help='Specify feature type in GTF annotation',
        default='gene',
        choices=['exon', 'gene'],
    )
    parser.add_argument('--genomeDir', help=HELP_DICT['genomeDir'])
    parser.add_argument('--featureCounts_param', help=HELP_DICT['additional_param'], default="")

    if sub_program:
        parser.add_argument('--input', help='Required. BAM file path.', required=True)
        parser = s_common(parser)
    return parser
