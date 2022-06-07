import ast
import argparse
import glob
import os

import pandas as pd
from plotnine import aes, geom_point, ggplot

from celescope.celescope import ArgFormatter
from celescope.tools import utils


SAMPLE_COL_INDEX = 2
MATCH_DIR_COL_INDEX = 3
VID_COL_INDEX = 4


@utils.add_log
def parse_mapfile(mapfile):
    sample_vid_dict = {}
    sample_match_dir_dict = {}
    df_mapfile = pd.read_csv(mapfile, sep='\t', header=None)

    def read_row(row):
        sample = row[SAMPLE_COL_INDEX]
        match_dir = row[MATCH_DIR_COL_INDEX]
        vid_list = [int(vid) for vid in str(row[VID_COL_INDEX]).strip().split(',')]
        sample_vid_dict[sample] = vid_list
        sample_match_dir_dict[sample] = match_dir

    df_mapfile.apply(read_row, axis=1)
    return sample_vid_dict, sample_match_dir_dict


class Plot_vid():
    def __init__(self, sample, outdir, vid_list, snp_dir, match_dir):
        self.sample = sample
        self.vid_list = vid_list

        # set
        vid_tsne_file = glob.glob(f'{snp_dir}/08.analysis_snp/*count_tsne.tsv')[0]
        self.df_vid_tsne = pd.read_csv(vid_tsne_file, sep='\t', converters={"VID": ast.literal_eval})
        match_tsne_file = glob.glob(f'{match_dir}/*analysis/*tsne_coord.tsv')[0]
        self.df_match_tsne = pd.read_csv(match_tsne_file, sep='\t', index_col=0)

        # out
        if not os.path.exists(outdir):
            os.system(f'mkdir -p {outdir}')
        self.out_prefix = f'{outdir}/{sample}'
        self.out_plot_file = f'{self.out_prefix}_VID_tsne.png'

    @utils.add_log
    def plot_vid(self):
        def set_label(row):
            for vid in self.vid_list:
                row["VIDs"] = "wild_type"
                if vid in row["VID"]:
                    row["VIDs"] = "mutation"
                    break
            return row
        df = self.df_vid_tsne.apply(set_label, axis=1)
        barcode_list = df.loc[df["VIDs"] == "mutation", ]["barcode"]
        self.df_match_tsne["VIDs"] = "wild_type"
        self.df_match_tsne.loc[barcode_list, "VIDs"] = "mutation"
        plot = ggplot(self.df_match_tsne, aes(x="tSNE_1", y="tSNE_2", color="VIDs")) + geom_point(size=0.2)
        plot.save(self.out_plot_file)


def main():
    parser = argparse.ArgumentParser(description='plot snp', formatter_class=ArgFormatter)
    parser.add_argument("--mapfile", help="mapfile with VIDs as 5th column", required=True)
    parser.add_argument("--outdir", help="output dir", default='plot_VID')
    args = parser.parse_args()

    sample_vid_dict, sample_match_dir_dict = parse_mapfile(args.mapfile)
    for sample in sample_vid_dict:
        vid_list = sample_vid_dict[sample]
        match_dir = sample_match_dir_dict[sample]

        runner = Plot_vid(
            sample=sample,
            outdir=args.outdir,
            vid_list=vid_list,
            snp_dir=sample,
            match_dir=match_dir
        )
        runner.plot_vid()


if __name__ == '__main__':
    main()
