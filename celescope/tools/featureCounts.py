
import os
import pathlib
from collections import defaultdict
from itertools import groupby
import subprocess
import unittest

import pandas as pd
import pysam

from celescope.rna.mkref import Mkref_rna
from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.tools import reference
from celescope.tools.__init__ import TAG_BAM_SUFFIX


GTF_TYPES = ['exon','gene']

def add_tag(seg, id_name, gene_correct_umi=None):
    """
    Args:
        seg: pysam bam segment
        id_name: {gene_id: gene_name}
        correct_dict: {low_seq: high_seq}

    Returns:
        seg with tag added

    Tags:
        CB: cell barcode
        UB: error-corrected UMI
        UR: original UMI
        GN: gene name
        GX: gene_id

    """
    attr = seg.query_name.split('_')
    barcode = attr[0]
    ur = ub = attr[1]

    # assign to some gene
    if seg.has_tag('XT'):
        gene_id = seg.get_tag('XT')
        # if multi-mapping reads are included in original bam,
        # there are multiple gene_ids
        if ',' in gene_id:
            gene_name = [id_name[i] for i in gene_id.split(',')]
            gene_name = ','.join(gene_name)
        else:
            gene_name = id_name[gene_id]
        seg.set_tag(tag='GN', value=gene_name, value_type='Z')
        seg.set_tag(tag='GX', value=gene_id, value_type='Z')

        if gene_correct_umi:
            if gene_id in gene_correct_umi and ur in gene_correct_umi[gene_id]:
                ub = gene_correct_umi[gene_id][ur]
    
    seg.set_tag(tag='CB', value=barcode, value_type='Z')
    seg.set_tag(tag='UB', value=ub, value_type='Z')
    seg.set_tag(tag='UR', value=ur, value_type='Z')

    return seg


def correct_umi(umi_dict, percent=0.1):
    """
    Correct umi_dict in place.
    Args:
        umi_dict: {umi_seq: umi_count}
        percent: if hamming_distance(low_seq, high_seq) == 1 and
            low_count / high_count < percent, merge low to high.
        return_dict: if set, return correct_dict = {low_seq:high_seq}
    Returns:
        n_corrected_umi: int
        n_corrected_read: int
    """
    n_corrected_umi = 0
    n_corrected_read = 0
    dic = {}

    # sort by value(UMI count) first, then key(UMI sequence)
    umi_arr = sorted(
        umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
    while True:
        # break when only highest in umi_arr
        if len(umi_arr) == 1:
            break
        umi_low = umi_arr.pop()
        low_seq = umi_low[0]
        low_count = umi_low[1]

        for umi_kv in umi_arr:
            high_seq = umi_kv[0]
            high_count = umi_kv[1]
            if float(low_count / high_count) > percent:
                break
            if utils.hamming_distance(low_seq, high_seq) == 1:
                n_low = umi_dict[low_seq]
                n_corrected_umi += 1
                n_corrected_read += n_low
                # merge
                umi_dict[high_seq] += n_low
                dic[low_seq] = high_seq
                del (umi_dict[low_seq])
                break

    return n_corrected_umi, n_corrected_read, dic


def discard_read(gene_umi_dict):
    """
    If two or more groups of reads have the same barcode and UMI, but different gene annotations, the gene annotation with the most supporting reads is kept for UMI counting, and the other read groups are discarded. In case of a tie for maximal read support, all read groups are discarded, as the gene cannot be confidently assigned.

    Returns:
        discarded_umi: set. umi with tie read count
        umi_gene_dict: {umi_seq: {gene_id: read_count}}
    """

    discard_umi = set()
    umi_gene_dict = defaultdict(lambda: defaultdict(int))
    for gene_id in gene_umi_dict:
        for umi in gene_umi_dict[gene_id]:
            umi_gene_dict[umi][gene_id] += gene_umi_dict[gene_id][umi]
    
    for umi in umi_gene_dict:
        max_read_count = max(umi_gene_dict[umi].values())
        gene_id_max = [gene_id for gene_id, read_count in umi_gene_dict[umi].items() if read_count==max_read_count]

        if len(gene_id_max) > 1:
            discard_umi.add(umi)
        else:
            gene_id = gene_id_max[0]
            umi_gene_dict[umi] = {gene_id: umi_gene_dict[umi][gene_id]}

    return discard_umi, umi_gene_dict


class FeatureCounts(Step):

    """
    ## Features
    - Assigning uniquely mapped reads to genomic features with FeatureCounts.
    ## Output
    - `{sample}` Numbers of reads assigned to features (or meta-features).
    - `{sample}_summary` Stat info for the overall summrization results, including number of 
    successfully assigned reads and number of reads that failed to be assigned due to 
    various reasons (these reasons are included in the stat info).
    - `{sample}_aligned_sortedByCoord_addTag.bam` featureCounts output BAM, 
    sorted by coordinates;BAM file contains tags as following(Software Version>=1.1.8):
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # set
        self.gtf = Mkref_rna.parse_genomeDir(self.args.genomeDir)['gtf']
        gp = reference.GtfParser(self.gtf)
        self.id_name = gp.get_id_name()

        # stats
        self.feature_log_dict = defaultdict(dict)
        self.n_corrected_read = 0
        self.n_corrected_umi = 0

        # out
        self.count_detail_file = f'{self.out_prefix}_count_detail.txt'
        input_basename = os.path.basename(self.args.input)
        self.featureCounts_bam = f'{self.outdir}/{input_basename}.featureCounts.bam'
        self.add_tag_bam = f'{self.out_prefix}_addTag.bam'
        self.out_bam = f'{self.out_prefix}_{TAG_BAM_SUFFIX}'


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
    def run_featureCounts(self, outdir, gtf_type):
        '''
        allow multimapping with -M; but each multi-mapped reads only have one alignment because of --outSAMmultNmax 1
        '''
        cmd = (
            'featureCounts '
            f'-s 1 '
            f'--largestOverlap '
            f'-M '
            f'-a {self.gtf} '
            f'-o {outdir}/{self.sample} '  
            '-R BAM '
            f'-T {self.args.thread} '
            f'-t {gtf_type} '
            f'{self.args.input} '
            '2>&1 '
        )
        if self.args.featureCounts_param:
            cmd += (" " + self.args.featureCounts_param)
        self.run_featureCounts.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run_get_log(self):
        tmp_dir = f'{self.outdir}/tmp/'
        for gtf_type in GTF_TYPES:
            outdir = f'{tmp_dir}/{gtf_type}'
            pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
            log_file = f'{outdir}/{self.sample}.summary'
            
            self.run_featureCounts(outdir, gtf_type)
            self.feature_log_dict[gtf_type] = FeatureCounts.read_log(log_file)

            if gtf_type == self.args.gtf_type:
                cmd = f'mv {outdir}/* {self.outdir}/ '
                subprocess.check_call(cmd, shell=True)

        cmd = f'rm -r {tmp_dir} '
        subprocess.check_call(cmd, shell=True)


    def remove_temp_file(self):
        os.remove(self.featureCounts_bam)
        os.remove(self.add_tag_bam)

    def run(self):
        self.run_get_log()
        self.add_metrics()
        self.get_count_detail_add_tag()
        utils.sort_bam(
            input_bam=self.add_tag_bam,
            output_bam=self.out_bam)
        self.remove_temp_file()

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

    @utils.add_log
    def get_count_detail_add_tag(self):
        """
        bam to detail table
        must be used on name_sorted bam
        Output file:
            - count_detail_file
            - bam with tag(remain name sorted)
        """
        save = pysam.set_verbosity(0)
        inputFile = pysam.AlignmentFile(self.featureCounts_bam, "rb")
        outputFile = pysam.AlignmentFile(self.add_tag_bam, 'wb', header=inputFile.header)
        pysam.set_verbosity(save)

        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'read', 'duplicate']) + '\n')

            def keyfunc(x):
                return x.query_name.split('_', 1)[0]
            for _, g in groupby(inputFile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                gene_umi_pos_cigar = utils.genDict(dim=4, valType=int)
                segs = []
                for seg in g:
                    segs.append(seg)
                    (barcode, umi) = seg.query_name.split('_')[:2]
                    if not seg.has_tag('XT'):
                        continue
                    gene_id = seg.get_tag('XT')
                    gene_umi_dict[gene_id][umi] += 1
                    gene_umi_pos_cigar[gene_id][umi][seg.reference_start][seg.cigarstring] += 1

                gene_correct_umi = None
                if self.args.correct_UMI:
                    gene_correct_umi = dict()
                    for gene_id in gene_umi_dict:
                        n_corrected_umi, n_corrected_read, correct_dict = correct_umi(gene_umi_dict[gene_id])
                        gene_correct_umi[gene_id] = correct_dict
                        self.n_corrected_read += n_corrected_read
                        self.n_corrected_umi += n_corrected_umi

                        # also correct umi in gene_umi_pos_cigar
                        for low_seq, high_seq in correct_dict.items():
                            for ref_start in gene_umi_pos_cigar[gene_id][low_seq]:
                                for cigar in gene_umi_pos_cigar[gene_id][low_seq][ref_start]:
                                    gene_umi_pos_cigar[gene_id][high_seq][ref_start][cigar] += 1
                            del gene_umi_pos_cigar[gene_id][low_seq]                    

                # output
                for gene_id in gene_umi_dict:
                    n_umi = len(gene_umi_dict[gene_id])
                    dup_list = []
                    n_read = 0
                    for umi in gene_umi_dict[gene_id]:
                        n_read += gene_umi_dict[gene_id][umi]
                        for pos in gene_umi_pos_cigar[gene_id][umi]:
                            for cigar in gene_umi_pos_cigar[gene_id][umi][pos]:
                                dup_list.append(str(gene_umi_pos_cigar[gene_id][umi][pos][cigar]))
                    duplicate = ','.join(dup_list)
                    fh1.write(f'{barcode}\t{gene_id}\t{n_umi}\t{n_read}\t{duplicate}\n')

                for seg in segs:
                    outputFile.write(add_tag(seg, self.id_name, gene_correct_umi))

        self.add_metric('n_corrected_read', self.n_corrected_read, show=False)
        self.add_metric('n_corrected_umi', self.n_corrected_umi, show=False)

        inputFile.close()
        outputFile.close()


@utils.add_log
def featureCounts(args):
    with FeatureCounts(args) as runner:
        runner.run()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument(
        '--gtf_type',
        help='Specify feature type in GTF annotation',
        default='gene',
        choices=['exon', 'gene'],
    )
    parser.add_argument('--genomeDir', help=HELP_DICT['genomeDir'])
    parser.add_argument('--featureCounts_param', help=HELP_DICT['additional_param'], default="")
    parser.add_argument('--correct_UMI', help='perform UMI correction.')

    if sub_program:
        parser.add_argument('--input', help='Required. BAM file path.', required=True)
        parser = s_common(parser)
    return parser


class featureCounts_test(unittest.TestCase):
    def test_correct_umi(self):
        dic = {
            "apple1": 2,
            "apple2": 30,
            "bears1": 5,
            "bears2": 10,
            "bears3": 100,
            "ccccc1": 20,
            "ccccc2": 199,
        }
        n_corrected_umi, n_corrected_read, _ = correct_umi(dic)
        dic_after_correct = {
            'ccccc1': 20,
            'apple2': 32,
            'bears3': 115,
            'ccccc2': 199,
        }
        self.assertEqual(dic, dic_after_correct)
        self.assertEqual(n_corrected_umi, 3)
        self.assertEqual(n_corrected_read, 2 + 5 + 10)

if __name__ == "__main__":
    unittest.main()