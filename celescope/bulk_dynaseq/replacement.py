"""Replacement step for bulk_dynaseq.

Processes each well's PosTag.bam independently (per-well background snp),
then merges all wells' results into a single h5ad/labeled/unlabeled output.
"""

import os
import sys
import pandas as pd
import pysam
from multiprocessing import Pool

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.dynaseq.__init__ import DYNA_MATRIX_DIR_SUFFIX
from celescope.tools.plotly_plot import Shared_box_plot
from celescope.tools.matrix import CountMatrix, Features


class Replacement(Step):
    """
    Features
    - Quantify unlabeled and labeled RNA per well, then merge into one output.
    - Boxplots for TOR rates distribution.

    Output
    - `{sample}.labeled.h5ad` h5ad file with ['total', 'labeled', 'unlabeled'] layers and TOR rate.
    - `labeled/` Merged labeled expression matrix.
    - `unlabeled/` Merged unlabeled expression matrix.
    - `{sample}_labeled_detail.csv` Detail file.
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # input files
        self.outdir = args.outdir
        self.sample = args.sample
        self.thread = int(args.thread)
        self.bamdir = args.bam_dir
        self.well_sample = utils.two_col_to_dict(args.well_sample)
        self.control_well = args.control_well
        self.bg = args.bg

        # set
        self.bg_external = []
        self.totaldf = pd.DataFrame()
        self.newdf = pd.DataFrame()
        self.olddf = pd.DataFrame()
        self.df = pd.DataFrame()

        # output files
        self.detail_txt = f"{self.out_prefix}_labeled_detail.csv"
        self.labeled_matrix = f"{self.out_prefix}_labeled_matrix.tsv.gz"
        self.unlabeled_matrix = f"{self.out_prefix}_unlabeled_matrix.tsv.gz"
        self.outs = [self.labeled_matrix, self.unlabeled_matrix]

    @utils.add_log
    def run(self):
        # get control well
        if self.bg is not None:
            self.bg_external.extend(self.bg)
        if self.control_well is not None:
            control_wells = self.control_well.split(',')
            control_samples = [self.well_sample[well] for well in control_wells]
            for sample in control_samples:
                if not os.path.exists(f"{self.bamdir}/{sample}/{sample}.snp.csv"):
                    print(f"Control sample {sample} has no snp file, skipped.")
                    continue
                self.bg_external.append(f"{self.bamdir}/{sample}/{sample}.snp.csv")
        else:
            control_samples = []

        # check bam existence, skip wells with no data
        valid_samples = []
        bam_list = []
        snp_list = []
        for well, sample in self.well_sample.items():
            if sample in control_samples:
                continue
            xbam = f"{self.bamdir}/{sample}/{sample}.PosTag.bam"
            xsnp = f"{self.bamdir}/{sample}/{sample}.snp.csv"
            xsnp_list = self.bg_external.copy()
            if not os.path.exists(xbam) or not os.path.exists(xsnp):
                print(f"Well {well}: {sample} has no data, skipped.")
                continue
            xsnp_list.append(xsnp)
            valid_samples.append(sample)
            bam_list.append(xbam)
            snp_list.append(xsnp_list)

        if len(bam_list) == 0:
            raise ValueError("No valid bam files found for replacement.")

        ## labeled and unlabeled
        self.dfs = self.run_quant(bam_list, snp_list, valid_samples)
        self.split_matrix()
        self.label_rate()
        self.quant_plot()
        self.quant_table()
        self.add_help()

    @utils.add_log
    def run_quant(self, bam_list, snp_list, sample_list):
        mincpu = min(len(sample_list), self.thread)
        with Pool(mincpu) as pool:
            results = pool.starmap(
                Replacement.quant, zip(bam_list, snp_list, sample_list)
            )
        return results

    def split_matrix(self):
        # merge matrix
        for i in self.dfs:
            self.totaldf = pd.concat([self.totaldf, i])
        self.newdf = self.totaldf[self.totaldf['TC'] == 1]
        self.olddf = self.totaldf[self.totaldf['TC'] == 0]
        outnewdf = self.newdf.groupby(['Well', 'gene']).agg({'UMI': 'count'})
        outolddf = self.olddf.groupby(['Well', 'gene']).agg({'UMI': 'count'})

        # output as matrix: rows=gene, columns=Well
        matrix_new = outnewdf['UMI'].unstack('Well').fillna(0).astype(int)
        matrix_old = outolddf['UMI'].unstack('Well').fillna(0).astype(int)
        matrix_new.index.name = ''
        matrix_old.index.name = ''
        matrix_new.to_csv(self.labeled_matrix, sep="\t")
        matrix_old.to_csv(self.unlabeled_matrix, sep="\t")

    @staticmethod
    def quant(bam_file, snp_list, sample):
        # get backgroud snp
        bg = Replacement.background_snp(snp_list)
        # get reads with TC
        df_total = Replacement.extract_dem(bam_file, bg, sample)
        return df_total

    @staticmethod
    def extract_dem(bam, bg, sample):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, 'rb')
        pysam.set_verbosity(save)
        countarr = []
        for read in bamfile.fetch():
            try:
                chro = read.reference_name
                ub = read.get_tag('UB')
                gene = read.get_tag('GX')

                if read.get_tag('ST') == '+':
                    stag = read.get_tag('TL')
                else:
                    stag = read.get_tag('AL')
                if len(stag) == 1 and stag[0] == 0:
                    tctag = 0
                else:
                    fcount = 0
                    for si in range(0, len(stag)):
                        pos = chro + '_' + str(stag[si])
                        if pos in bg:
                            fcount += 1
                    if fcount == len(stag):
                        tctag = 0
                    else:
                        tctag = 1
                countarr.append([gene, sample, ub, tctag])
            except (ValueError, KeyError):
                continue
        bamfile.close()

        totaldf = pd.DataFrame(countarr)
        totaldf.columns = ['gene', 'Well', 'UMI', 'TC']
        df_unique = totaldf.groupby(['gene', 'Well', 'UMI']).apply(
            lambda x: x.loc[x['TC'].idxmax()]
        ).reset_index(drop=True)

        return df_unique

    @staticmethod
    def background_snp(bgfiles):
        outdict = {}
        for bgfile in bgfiles:
            if bgfile.endswith('.csv'):
                df = pd.read_csv(bgfile, dtype={"chrom": str})
                if df.shape[0] == 0:
                    continue
                df['chrpos'] = df['chrom'] + '_' + df['pos'].astype(str)
                df1 = df[['chrpos', 'posratio']]
                df1.set_index('chrpos', inplace=True)
                if len(df1.index) > 0:
                    outdict.update(df1.to_dict(orient='index'))
            elif bgfile.endswith('.vcf'):
                bcf_in = pysam.VariantFile(bgfile)
                for rec in bcf_in.fetch():
                    try:
                        chrom, pos = rec.chrom, rec.pos
                        chr_pos = chrom + '_' + str(pos - 1)
                        outdict[chr_pos] = 1
                    except (ValueError, KeyError):
                        continue
                bcf_in.close()
            else:
                sys.exit(
                    'Background snp file format cannot be recognized! Only csv or vcf format.'
                )
        return outdict

    @utils.add_log
    def label_rate(self):
        df_sum = self.totaldf.groupby('Well').agg({
            'UMI': 'count',
            'gene': 'nunique'
        })
        df_new_sum = self.newdf.groupby('Well').agg({
            'UMI': 'count',
            'gene': 'nunique'
        })
        self.df = pd.concat([df_sum, df_new_sum], axis=1)
        self.df.columns = ['UMI', 'gene', 'labeled_UMI', 'labeled_gene']
        self.df.loc[:, 'UMI_rate'] = self.df['labeled_UMI'] / self.df['UMI']
        self.df.loc[:, 'gene_rate'] = self.df['labeled_gene'] / self.df['gene']
        self.df = self.df[['UMI_rate', 'gene_rate']]
        self.df.to_csv(self.detail_txt, sep='\t')

    @utils.add_log
    def quant_plot(self):
        outplot = self.df.copy()
        quantplot = Shared_box_plot(df_box=outplot).get_plotly_div()
        self.add_data(quant_box=quantplot)

    @utils.add_log
    def quant_table(self):
        table_dict = self.get_table_dict(
            title='',
            table_id='bulk_replacement',
            df_table=self.df.reset_index(),
        )
        self.add_data(table_dict=table_dict)

    @utils.add_log
    def add_help(self):
        self.add_help_content(
            name='labeled_rate',
            content='Fraction of labeled UMI rate and gene label rate per well'
        )


@utils.add_log
def replacement(args):
    with Replacement(args, display_title="Labeled") as runner:
        runner.run()


def get_opts_replacement(parser, sub_program):
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"])
    parser.add_argument(
        "--control_well",
        help="Control well number, split by comma. Default is None",
        required=False,
        type=str,
        default=None,
    )
    parser.add_argument(
        "--replacementMem",
        default=10,
        type=int,
        help="Set replacement memory.",
    )
    parser.add_argument(
        "--well_sample",
        help="tsv file of well numbers and sample names. The first column is well numbers and the second column is sample names.",
        required=True,
    )
    if sub_program:
        parser.add_argument(
            "--bam_dir",
            help="Directory containing per-well conversion outputs",
            required=True,
        )
        parser.add_argument(
            "--bg",
            nargs="+",
            required=False,
            help="external background snp file, csv or vcf format",
        )
        parser = s_common(parser)
    return parser
