import os

import pandas as pd
import pysam

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH, HELP_DICT
from celescope.fusion.mkref import parse_genomeDir_fusion
from celescope.tools.step import Step, s_common



def is_fusion(pos, read_start, read_length, flanking_base):
    test_start = (pos - flanking_base) >= read_start
    test_end = (pos + flanking_base) <= (read_start + read_length)
    return test_start and test_end


class CountFusion(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.flanking_base = int(args.flanking_base)
        self.UMI_min = int(args.UMI_min)
        self.match_dir = args.match_dir
        self.fusion_genomeDir = args.fusion_genomeDir
        self.fusion_pos_file = parse_genomeDir_fusion(self.fusion_genomeDir)['fusion_pos']
        self.bam = args.bam

        # out
        self.out_read_count_file = self.out_prefix + "_fusion_read_count.tsv"
        self.out_umi_count_file = self.out_prefix + "_fusion_UMI_count.tsv"
        self.out_barcode_count_file = self.out_prefix + "_fusion_barcode_count.tsv"
        self.out_tsne_file = self.out_prefix + "_fusion_tsne.tsv"

    def read_pos(self):
        dic = {}
        df = pd.read_csv(self.fusion_pos_file, sep="\t")
        for tag, pos in zip(df.iloc[:, 0], df.iloc[:, 1]):
            dic[tag] = pos
        return dic

    def count_fusion(self):
        fusion_pos = self.read_pos()
        # barcode
        match_barcode, _n_barcode = utils.read_barcode_file(self.match_dir)
        # tsne
        match_tsne_file = utils.parse_match_dir(self.match_dir)["tsne_coord"]
        df_tsne = pd.read_csv(match_tsne_file, sep="\t", index_col=0)

        # process bam
        samfile = pysam.AlignmentFile(self.bam, "rb")
        header = samfile.header
        new_bam = pysam.AlignmentFile(self.out_prefix + "_fusion.bam", "wb", header=header)
        count_dic = utils.genDict(dim=3)
        for read in samfile:
            tag = read.reference_name
            read_start = int(read.reference_start)
            read_length = len(read.query_sequence)
            attr = read.query_name.split("_")
            barcode = attr[0]
            umi = attr[1]
            if tag in fusion_pos.keys():
                if barcode in match_barcode:
                    if is_fusion(
                        pos=fusion_pos[tag],
                        read_start=read_start,
                        read_length=read_length,
                        flanking_base=self.flanking_base,
                    ):
                        new_bam.write(read)
                        count_dic[barcode][tag][umi] += 1
        new_bam.close()

        # write dic to pandas df
        rows = []
        for barcode in count_dic:
            for tag in count_dic[barcode]:
                for umi in count_dic[barcode][tag]:
                    rows.append([barcode, tag, umi, count_dic[barcode][tag][umi]])
        df_read = pd.DataFrame(rows)
        df_read.rename(
            columns={0: "barcode", 1: "tag", 2: "UMI", 3: "read_count"}, inplace=True
        )
        df_read.to_csv(self.out_read_count_file, sep="\t", index=False)

        if not rows:
            count_fusion.logger.error("***** NO FUSION FOUND! *****")
        else:
            df_umi = df_read.groupby(["barcode", "tag"]).agg({"UMI": "count"})
            df_umi = df_umi[df_umi["UMI"] >= self.UMI_min]
            df_umi.to_csv(self.out_umi_count_file, sep="\t")

            df_umi.reset_index(inplace=True)
            df_barcode = df_umi.groupby(["tag"]).agg({"barcode": "count"})
            n_match_barcode = len(match_barcode)
            # add zero count tag
            for tag in fusion_pos.keys():
                if not tag in df_barcode.barcode:
                    new_row = pd.Series(data={"barcode": 0}, name=tag)
                    df_barcode = df_barcode.append(new_row, ignore_index=False)
            df_barcode["percent"] = df_barcode["barcode"] / n_match_barcode
            df_barcode.to_csv(self.out_barcode_count_file, sep="\t")

            df_pivot = df_umi.pivot(index="barcode", columns="tag", values="UMI")
            df_pivot.fillna(0, inplace=True)
            df_tsne_fusion = pd.merge(
                df_tsne, df_pivot, right_index=True, left_index=True, how="left"
            )
            df_tsne_fusion.fillna(0, inplace=True)
            df_tsne_fusion.to_csv(self.out_tsne_file, sep="\t")

            # plot
            count_fusion.logger.info("plot fusion...!")
            app = f'{ROOT_PATH}/fusion/plot_fusion.R'
            cmd = f"Rscript {app} --tsne_fusion {self.out_tsne_file} --outdir {self.outdir}"
            os.system(cmd)
            count_fusion.logger.info("plot done.")

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
