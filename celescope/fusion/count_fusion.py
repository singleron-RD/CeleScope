import subprocess

import pandas as pd
import pysam

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH, HELP_DICT
from celescope.fusion.mkref import parse_genomeDir_fusion
from celescope.tools.step import Step, s_common


class CountFusion(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.flanking_base = int(args.flanking_base)
        self.UMI_min = int(args.UMI_min)
        self.match_dir = args.match_dir
        self.fusion_genomeDir = args.fusion_genomeDir
        self.fusion_pos_file = parse_genomeDir_fusion(self.fusion_genomeDir)['fusion_pos']
        self.bam = args.bam

        # set
        self.pos_dic = self.read_pos()
        self.match_barcode, self.n_match_barcode = utils.read_barcode_file(self.match_dir)
        # tsne
        match_tsne_file = utils.parse_match_dir(self.match_dir)["tsne_coord"]
        self.df_tsne = pd.read_csv(match_tsne_file, sep="\t", index_col=0)

        # out
        self.fusion_bam = self.out_prefix + "_fusion.bam"
        self.out_read_count_file = self.out_prefix + "_fusion_read_count.tsv"
        self.out_umi_count_file = self.out_prefix + "_fusion_UMI_count.tsv"
        self.out_barcode_count_file = self.out_prefix + "_fusion_barcode_count.tsv"
        self.out_tsne_file = self.out_prefix + "_fusion_tsne.tsv"

    def read_pos(self):
        """
        pos_dic
            key: sequence name in fusion.fasta
            value: end postion of the first gene(1-based).
        """
        pos_dic = {}
        df = pd.read_csv(self.fusion_pos_file, sep="\t")
        for name, pos in zip(df.iloc[:, 0], df.iloc[:, 1]):
            pos_dic[name] = pos
        return pos_dic

    def count_fusion(self):
        """
        find valid fusion reads
            1. flank the fusion position
            2. match barcode
        """
        count_dic = utils.genDict(dim=3)

        with pysam.AlignmentFile(self.bam, "rb") as bam:
            header = bam.header
            with pysam.AlignmentFile(self.fusion_bam, "wb", header=header) as fusion_bam:
                for name in self.pos_dic:
                    pos = self.pos_dic[name]
                    left = pos - self.flanking_base
                    right = pos + self.flanking_base
                    for read in bam.fetch(
                        reference=name,
                        start=left,
                        end=right,
                    ):
                        left_bases = read.get_overlap(left, pos)
                        right_bases = read.get_overlap(pos, right)
                        if left_bases < self.flanking_base or right_bases <  self.flanking_base:
                            continue
                        attr = read.query_name.split("_")
                        barcode = attr[0]
                        umi = attr[1]
                        if barcode in self.match_barcode:
                            fusion_bam.write(read)
                            count_dic[barcode][name][umi] += 1

        # write dic to pandas df
        rows = []
        for barcode in count_dic:
            for name in count_dic[barcode]:
                for umi in count_dic[barcode][name]:
                    rows.append([barcode, name, umi, count_dic[barcode][name][umi]])
        if not rows:
            count_fusion.logger.warning("***** NO FUSION FOUND! *****")
            return
        df_read = pd.DataFrame(rows)
        df_read.rename(
            columns={0: "barcode", 1: "name", 2: "UMI", 3: "read_count"}, inplace=True
        )
        df_read.to_csv(self.out_read_count_file, sep="\t", index=False)

        df_umi = df_read.groupby(["barcode", "name"]).agg({"UMI": "count"})
        df_umi = df_umi[df_umi["UMI"] >= self.UMI_min]
        df_umi.to_csv(self.out_umi_count_file, sep="\t")

        df_umi.reset_index(inplace=True)
        df_barcode = df_umi.groupby(["name"]).agg({"barcode": "count"})
        # add zero count name
        for name in self.pos_dic:
            if not name in df_barcode.index:
                new_row = pd.Series(data={"barcode": 0}, name=name)
                df_barcode = df_barcode.append(new_row, ignore_index=False)
        df_barcode["percent"] = df_barcode["barcode"] / self.n_match_barcode
        df_barcode.to_csv(self.out_barcode_count_file, sep="\t")

        df_pivot = df_umi.pivot(index="barcode", columns="name", values="UMI")
        df_pivot.fillna(0, inplace=True)
        df_tsne_fusion = pd.merge(
            self.df_tsne, df_pivot, right_index=True, left_index=True, how="left"
        )
        df_tsne_fusion.fillna(0, inplace=True)
        df_tsne_fusion.to_csv(self.out_tsne_file, sep="\t")

        self.plot_fusion()

    @utils.add_log
    def plot_fusion(self):
        app = f'{ROOT_PATH}/fusion/plot_fusion.R'
        cmd = f"Rscript {app} --tsne_fusion {self.out_tsne_file} --outdir {self.outdir}"
        subprocess.check_call(cmd, shell=True)

    def run(self):
        self.count_fusion()
        self.clean_up()


@utils.add_log
def count_fusion(args):
    step_name = "count_fusion"
    runner = CountFusion(args, step_name)
    runner.run()


def get_opts_count_fusion(parser, sub_program):
    if sub_program:
        s_common(parser)
        parser.add_argument("--bam", help='STAR bam file.', required=True)
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'], required=True)
    parser.add_argument('--fusion_genomeDir', help='Fusion genome directory.', required=True)
    parser.add_argument(
        "--flanking_base", 
        help="Number of bases flanking the fusion position.",
        default=5)
    parser.add_argument(
        "--UMI_min",
        help="Minimum number of fusion UMI to consider a cell as a cell with fusion event.",
        default=1)
