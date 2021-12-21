import argparse
import glob
import os
import subprocess

import pandas as pd
from plotnine import ggplot, aes, geom_line

from celescope.celescope import ArgFormatter
from celescope.__init__ import HELP_DICT, ROOT_PATH
from celescope.rna.mkref import parse_genomeDir_rna
import celescope.tools.utils as utils

SAMPLE_COL_INDEX = 2


def parse_mapfile(mapfile):
    sample_set = set()
    df_mapfile = pd.read_csv(mapfile, sep='\t', header=None)

    def read_row(row):
        sample = row[SAMPLE_COL_INDEX]
        sample_set.add(sample)

    df_mapfile.apply(read_row, axis=1)
    return sample_set


class Mt_summary():
    def __init__(self, sample, outdir, genomeDir, root_dir):
        self.sample = sample
        self.outdir = outdir

        # set
        match_dir = f'{root_dir}/{sample}'
        self.mt_gene_list_file = parse_genomeDir_rna(genomeDir)['mt_gene_list']
        _barcodes, self.ncell = utils.read_barcode_file(match_dir)
        self.bam = None
        try:
            self.bam = glob.glob(
                f'{match_dir}/03*/{sample}*sortedByCoord.out.bam')[0]
        except IndexError:
            print("STAR bam does not exist! Skip coverage summary.")

        self.matrix_dir = glob.glob(f'{match_dir}/*count/{sample}_matrix_10X')[0]

        # out
        if not os.path.exists(outdir):
            os.system(f'mkdir -p {outdir}')
        out_prefix = f'{outdir}/{sample}'
        self.mt_bam = f'{out_prefix}_mt.bam'
        self.mt_depth = f'{out_prefix}_mt_depth.tsv'
        self.coverage_plot = f'{out_prefix}_mt_coverage.png'

    @utils.add_log
    def samtools(self):
        cmd = (
            f'samtools view -b {self.bam} MT -o {self.mt_bam};'
            f'samtools depth -a {self.mt_bam} > {self.mt_depth}'
        )
        self.samtools.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def umi_summary(self):
        cmd = (
            f'Rscript {ROOT_PATH}/scripts/gene_umi_summary.R '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
            f'--mt_gene_list_file {self.mt_gene_list_file} '
            f'--matrix_dir {self.matrix_dir} '
        )
        self.umi_summary.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def coverage_summary(self):
        self.samtools()
        df = pd.read_csv(self.mt_depth, sep='\t', header=None)
        df.columns = ["MT", "position", "read_count"]
        df["mean_read_count_per_cell"] = df["read_count"].apply(lambda x: x / self.ncell)
        plot = ggplot(df, aes(x="position", y="mean_read_count_per_cell")) + geom_line()
        plot.save(self.coverage_plot)

    @utils.add_log
    def run(self):
        self.umi_summary()
        if self.bam:
            self.coverage_summary()


def main():
    parser = argparse.ArgumentParser(description='plot snp', formatter_class=ArgFormatter)
    parser.add_argument("--mapfile", help="mapfile with VIDs as 5th column", required=True)
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"],
                        default='/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92')
    parser.add_argument("--root_dir", help='input root_dir', default='./')
    parser.add_argument("--outdir", help="output dir", default='mt_summary')
    args = parser.parse_args()

    sample_set = parse_mapfile(args.mapfile)
    for sample in sample_set:
        runner = Mt_summary(
            sample=sample,
            outdir=args.outdir,
            genomeDir=args.genomeDir,
            root_dir=args.root_dir,
        )
        runner.run()


if __name__ == '__main__':
    main()
