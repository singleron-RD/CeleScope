from collections import defaultdict
import itertools

import pysam
import pandas as pd
from celescope.tools.step import Step

from celescope.tools.count import Count
from celescope.tools import utils
from celescope.tools.step import Step, s_common


class Count_capture_virus_mtx(Count):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.bam = args.bam
        self.gtf = args.gtf
        self.filter_umi_file = args.filter_umi_file
        self.count_detail_file = f'{self.outdir}/{self.sample}_count_detail.txt'
        self.matrix_dir = f'{self.outdir}/{self.sample}_virus_matrix'
        self.id_name = utils.get_id_name_dict(self.gtf)

    @staticmethod
    def get_barcodes(filter_umi_file):
        df = pd.read_table(filter_umi_file)
        return set(list(df['barcode']))

    @staticmethod
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

    def bam2table(self):
        """
        bam to detail table
        must be used on name_sorted bam
        """
        samfile = pysam.AlignmentFile(self.bam, "rb")
        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'count']) + '\n')

            def keyfunc(x):
                return x.query_name.split('_', maxsplit=1)[0]
            for _, g in itertools.groupby(samfile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                for seg in g:
                    (barcode, umi) = seg.query_name.split('_')[:2]
                    if not seg.has_tag('XT'):
                        continue
                    gene_id = seg.get_tag('XT')
                    gene_umi_dict[gene_id][umi] += 1
                for gene_id in gene_umi_dict:
                    Count.correct_umi(gene_umi_dict[gene_id])

                discard_umi, umi_gene_dict = Count_capture_virus_mtx.discard_read(gene_umi_dict)

                # output
                for umi in umi_gene_dict:
                    if umi not in discard_umi:
                        for gene_id in umi_gene_dict[umi]:
                            fh1.write('%s\t%s\t%s\t%s\n' % (barcode, gene_id, umi,
                                                            umi_gene_dict[umi][gene_id]))
        samfile.close()
    
    
    

    def run(self):
        self.bam2table()
        df = pd.read_table(self.count_detail_file, header=0)
        barcodes = Count_capture_virus_mtx.get_barcodes(self.args.filter_umi_file)
        df_filter = df[df["Barcode"].isin(barcodes)]
        df_filter.index = df["Barcode"].to_list()
        self.write_matrix_10X(df_filter, self.matrix_dir)


@utils.add_log
def count_capture_virus_mtx(args):
    step_name = 'count_capture_virus_mtx'
    runner = Count_capture_virus_mtx(args, step_name)
    runner.run()

def get_opts_count_capture_virus_mtx(parser, sub_program):
    parser.add_argument(
        "--gtf",
        help="Optional. Genome gtf file. Use absolute path or relative path to `genomeDir`.",
        )
    if sub_program:
        parser.add_argument('--filter_umi_file', help='filter umi file', required=True)
        parser.add_argument('--bam', help='Required. BAM file from featureCounts.', required=True)
        s_common(parser)