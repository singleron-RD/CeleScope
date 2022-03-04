import pysam

import pandas as pd

from celescope.tools.capture.count_bam import Count_bam, get_opts_count_bam
from celescope.fusion.mkref import Mkref_fusion



class Count_fusion(Count_bam):
    def __init__(self, args, display_title='Count'):
        super().__init__(args, display_title)

        self.flanking_base = int(args.flanking_base)
        fusion_pos_file = Mkref_fusion.parse_genomeDir(args.fusion_genomeDir)['fusion_pos']
        self.pos_dict = self.read_pos_file(fusion_pos_file)
        self.fusion_bam = f'{self.out_prefix}_raw_fusion.bam'
 
    @staticmethod
    def read_pos_file(fusion_pos_file):
        """
        pos_dic
            key: sequence name in fusion.fasta
            value: end postion of the first gene(1-based).
        """
        pos_dict = {}
        df = pd.read_csv(fusion_pos_file, sep="\t")
        for name, pos in zip(df.iloc[:, 0], df.iloc[:, 1]):
            pos_dict[name] = pos
        return pos_dict

    def process_bam(self):
        """
        find valid fusion reads
            1. flank the fusion position
            2. match barcode
        """

        with pysam.AlignmentFile(self.capture_bam, "rb") as bam:
            header = bam.header
            with pysam.AlignmentFile(self.fusion_bam, "wb", header=header) as fusion_bam:
                for ref in self.pos_dict:
                    pos = self.pos_dict[ref]
                    left = pos - self.flanking_base
                    right = pos + self.flanking_base
                    for read in bam.fetch(
                        reference=ref,
                        start=left,
                        end=right,
                    ):
                        left_bases = read.get_overlap(left, pos)
                        right_bases = read.get_overlap(pos, right)
                        if left_bases < self.flanking_base or right_bases < self.flanking_base:
                            continue
                        attr = read.query_name.split("_")
                        barcode = attr[0]
                        umi = attr[1]
                        if barcode in self.match_barcode:
                            fusion_bam.write(read)
                            self.count_dict[barcode][ref][umi] += 1


def count_fusion(args):
    with Count_fusion(args) as runner:
        runner.run()


def get_opts_count_fusion(parser, sub_program):
    parser.add_argument('--fusion_genomeDir', help='Fusion genome directory.', required=True)
    parser.add_argument(
        "--flanking_base",
        help="Number of bases flanking the fusion position.",
        default=5)
    get_opts_count_bam(parser, sub_program)

