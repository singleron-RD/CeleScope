from collections import defaultdict
import itertools

import pysam
import pandas as pd

from celescope.tools.featureCounts import correct_umi, discard_read
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools import reference
from celescope.tools.matrix import CountMatrix


class Count(Step):
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)
        self.bam = args.bam
        gp = reference.GtfParser(args.gtf)
        gp.get_id_name()
        self.features = gp.get_features()

        self.count_detail_file = f"{self.outdir}/{self.sample}_count_detail.txt"
        self.matrix_dir = f"{self.outdir}/{self.sample}_virus_matrix"

    @staticmethod
    def get_barcodes(filter_umi_file):
        df = pd.read_csv(filter_umi_file)
        return list(df["barcode"])

    def run(self):
        self.bam2table()
        df = pd.read_table(self.count_detail_file, header=0)
        barcodes = Count.get_barcodes(self.args.filter_umi_file)
        df = df.groupby(["Barcode", "geneID"]).agg({"UMI": "count"})

        count_matrix = CountMatrix.from_dataframe(
            df, self.features, barcodes=barcodes, value="UMI"
        )
        count_matrix.to_matrix_dir(self.matrix_dir)

    @utils.add_log
    def bam2table(self):
        """
        bam to detail table
        must be used on name_sorted bam
        """
        samfile = pysam.AlignmentFile(self.bam, "rb")
        with open(self.count_detail_file, "wt") as fh1:
            fh1.write("\t".join(["Barcode", "geneID", "UMI", "count"]) + "\n")

            def keyfunc(x):
                return x.query_name.split(":", maxsplit=1)[0]

            for _, g in itertools.groupby(samfile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                for seg in g:
                    (barcode, umi) = seg.query_name.split(":")[:2]
                    if not seg.has_tag("XT"):
                        continue
                    gene_id = seg.get_tag("XT")
                    gene_umi_dict[gene_id][umi] += 1
                for gene_id in gene_umi_dict:
                    correct_umi(gene_umi_dict[gene_id])

                discard_umi, umi_gene_dict = discard_read(gene_umi_dict)

                # output
                for umi in umi_gene_dict:
                    if umi not in discard_umi:
                        for gene_id in umi_gene_dict[umi]:
                            fh1.write(
                                "%s\t%s\t%s\t%s\n"
                                % (barcode, gene_id, umi, umi_gene_dict[umi][gene_id])
                            )
        samfile.close()


@utils.add_log
def count(args):
    if args.gtf:
        with Count(args, display_title="Cells") as runner:
            runner.run()


def get_opts_count(parser, sub_program):
    parser.add_argument(
        "--gtf",
        help="Optional. Genome gtf file. Use absolute path or relative path to `genomeDir`.",
    )
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--filter_umi_file", help="filter umi file", required=True)
        parser.add_argument(
            "--bam", help="Required. BAM file from featureCounts.", required=True
        )
