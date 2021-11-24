
import os
import re
import subprocess

import pysam
import celescope.tools.utils as utils
from celescope.rna.mkref import parse_genomeDir_rna
from celescope.tools.step import Step, s_common


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

    def __init__(self, args,display_title=None):
        Step.__init__(self, args,display_title=display_title)

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
            self.add_metric(
                name='Assigned', 
                value=tmp_arr[0],
                help_info='reads that can be successfully assigned without ambiguity'
            )
            self.add_metric(
                name='Unassigned_NoFeatures', 
                value=tmp_arr[0],
                help_info='alignments that do not overlap any feature'
            )
            self.add_metric(
                name='Unassigned_Ambiguity', 
                value=tmp_arr[0],
                help_info='alignments that overlap two or more features'
            )         
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

    @utils.add_log
    def name_sort_bam(self):
        cmd = (
            'samtools sort -n '
            f'-@ {self.thread} '
            f'-o {self.name_sorted_bam} '
            f'{self.featureCounts_bam} '
        )
        FeatureCounts.name_sort_bam.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.run_featureCounts()
        add_tag(self.featureCounts_bam, self.gtf)
        self.name_sort_bam()
        self.format_stat() 


@utils.add_log
def add_tag(bam, gtf):
    """
    - CB cell barcode
    - UB UMI
    - GN gene name
    - GX gene id
    """
    id_name = utils.get_id_name_dict(gtf)
    temp_bam_file = bam + ".temp"

    with pysam.AlignmentFile(bam, "rb") as original_bam:
        header = original_bam.header
        with pysam.AlignmentFile(temp_bam_file, "wb", header=header) as temp_bam:
            for read in original_bam:
                attr = read.query_name.split('_')
                barcode = attr[0]
                umi = attr[1]
                read.set_tag(tag='CB', value=barcode, value_type='Z')
                read.set_tag(tag='UB', value=umi, value_type='Z')
                # assign to some gene
                if read.has_tag('XT'):
                    gene_id = read.get_tag('XT')
                    # if multi-mapping reads are included in original bam,
                    # there are multiple gene_ids
                    if ',' in gene_id:
                        gene_name = [id_name[i] for i in gene_id.split(',')]
                        gene_name = ','.join(gene_name)
                    else:
                        gene_name = id_name[gene_id]
                    read.set_tag(tag='GN', value=gene_name, value_type='Z')
                    read.set_tag(tag='GX', value=gene_id, value_type='Z')
                temp_bam.write(read)
    cmd = f'mv {temp_bam_file} {bam}'
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def featureCounts(args):
    with FeatureCounts(args) as runner:
        runner.run()


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
