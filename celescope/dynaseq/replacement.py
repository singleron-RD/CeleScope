import os
import re
import glob
import pandas as pd
import scanpy as sc
import pysam
from multiprocessing import Pool
import subprocess

from celescope.tools.step import Step, s_common
from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.dynaseq.__init__ import DYNA_MATRIX_DIR_SUFFIX
from celescope.tools.matrix import CountMatrix, Features, ROW, COLUMN
from celescope.tools.plotly_plot import Tsne_plot, Violin_plot

toolsdir = os.path.dirname(__file__)


class Replacement(Step):
    """
    Features
    - Quantify unlabeled and labeled RNA.
    - Boxplots for TOR rates distribution.
    - TSNE plot for TOR rate

    Output
    - `{sample}.labeled.h5ad` h5ad file contains ['total', 'labeled', 'unlabeled'] layers and TOR rate of each cell/gene.
    - `{sample}_labeled_feature_bc_matrix` The labeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
    - `{sample}_unlabeled_feature_bc_matrix` The unlabeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
    - `{sample}_labeled_detail.txt`  tab-delimited  file:
        - Barcode: Cell barcode sequence
        - UMI: UMI sequence
        - geneID: gene ID
        - TC: TC site number in a read (backgroup snp removed)
    """

    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        # input files
        self.outdir = args.outdir
        self.sample = args.sample
        self.thread = int(args.thread)
        self.bcfile = args.cell
        self.featurefile = args.gene
        self.snp_file = args.bg
        self.tsne = args.tsne

        # set
        if args.bam:
            self.inbam = args.bam
        else:
            self.inbam = None
        self.cell_list, self.cell_num = utils.read_one_col(self.bcfile)
        self.features = Features.from_tsv(tsv_file=self.featurefile)
        self.totaldf = pd.DataFrame()
        self.newdf, self.olddf = pd.DataFrame(), pd.DataFrame()
        matrix_file = os.path.dirname(self.bcfile)
        self.adata = sc.read_10x_mtx(
            matrix_file,
            var_names="gene_ids",
        )
        self.bg = None

        # output files
        ## tmp outputdir
        self.tmp_dir = f"{args.outdir}/../tmp/"
        utils.check_mkdir(self.tmp_dir)
        ## final outputs
        self.h5ad = f"{self.out_prefix}.labeled.h5ad"
        self.detail_txt = f"{self.out_prefix}_labeled_detail.csv"
        self.dir_labeled = f"{self.outdir}/{DYNA_MATRIX_DIR_SUFFIX[0]}"
        self.dir_unlabeled = f"{self.outdir}/{DYNA_MATRIX_DIR_SUFFIX[1]}"
        self.outs = [self.dir_labeled, self.dir_unlabeled, self.h5ad]

    @utils.add_log
    def run(self):
        # get backgroud snp
        self.bg = self.background_snp()
        # replacement
        dfs = self.run_quant()
        self.labeled(dfs)
        # stat and plot
        self.tor_plot()
        self.add_help()
        # output dedup and clean
        self.clean_tmp()

    @utils.add_log
    def run_quant(self):
        # tmp_dir = f"{self.tmp_dir}/2"
        # utils.check_mkdir(tmp_dir)
        ## set Parallelism para
        fetch_arr = glob.glob(f"{self.tmp_dir}/1/*.bam", recursive=False)
        if len(fetch_arr) == 0:
            if self.inbam:
                fetch_arr = [self.inbam]
            else:
                raise ValueError(
                    "No bam detected, please check `../tmp/1` or specific by `--bam`."
                )
        snp_list = [self.bg] * len(fetch_arr)
        cell_arr = [self.cell_list] * len(fetch_arr)

        mincpu = min(self.cell_num, self.thread)
        with Pool(mincpu) as pool:
            results = pool.starmap(
                Replacement.modify_bam, zip(fetch_arr, snp_list, cell_arr)
            )
        return results

    @utils.add_log
    def labeled(self, dfs):
        # merge df
        df = pd.concat(dfs, ignore_index=True)
        self.totaldf = df.loc[df.groupby(["Barcode", "geneID", "UMI"])["TC"].idxmax()]
        self.newdf = self.totaldf[self.totaldf["TC"] > 0]
        self.olddf = self.totaldf[self.totaldf["TC"] == 0]
        self.totaldf.to_csv(self.detail_txt, index=False)

        # output
        tmp_newdf = self.newdf.groupby([COLUMN, ROW]).agg({"UMI": "count"})
        tmp_olddf = self.olddf.groupby([COLUMN, ROW]).agg({"UMI": "count"})
        self.write_sparse_matrix(tmp_newdf, self.dir_labeled)
        self.write_sparse_matrix(tmp_olddf, self.dir_unlabeled)
        self.write_h5ad()

    @staticmethod
    def modify_bam(bam, bg, cells):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, "rb")
        pysam.set_verbosity(save)
        readdict = {}

        for read in bamfile.fetch(until_eof=True):
            try:
                cb = read.get_tag("CB")
                if cb not in cells:
                    continue
                chro = read.reference_name
                ub = read.get_tag("UB")
                gene = read.get_tag("GX")
                tctag = 0
                true_tc = []

                if read.get_tag("ST") == "+":
                    stag = read.get_tag("TL")
                else:
                    stag = read.get_tag("AL")
                if len(stag) == 1 and stag[0] == 0:
                    tctag = 0
                    true_tc = stag
                else:
                    for si in range(0, len(stag)):
                        pos = chro + "_" + str(stag[si])
                        if pos not in bg:
                            true_tc.append(int(stag[si]))
                    tctag = len(true_tc)
                ## dedup: select the most TC read per UMI_gene
                readid = ":".join([cb, ub, gene])
                if readid not in readdict:
                    readdict[readid] = tctag
                else:
                    if tctag > readdict[readid]:
                        readdict[readid] = tctag

            except (ValueError, KeyError):
                continue
        bamfile.close()

        ## count df
        if len(readdict) == 0:
            return pd.DataFrame()
        tc_df = pd.DataFrame.from_dict(readdict, orient="index", columns=["TC"])
        tc_df["TC"] = tc_df["TC"].astype("uint8")
        tc_df.reset_index(inplace=True)
        ub_df = tc_df["index"].str.split(":", expand=True)
        ub_df.columns = ["Barcode", "UMI", "geneID"]
        tc_df = pd.concat([ub_df, tc_df["TC"]], axis=1)
        return tc_df

    @staticmethod
    def createTag(d):
        return "".join(["".join(key) + str(d[key]) + ";" for key in d.keys()])[:-1]

    @staticmethod
    def modifySCTag(sc, cnt, sstag):
        pattern = re.compile(rf"{sstag}\d+")
        result = pattern.sub(f"{sstag}{cnt}", sc)
        return result

    @utils.add_log
    def background_snp(self):
        outdict = {}
        bgs = []
        for bgargv in self.snp_file:
            if "," in bgargv:
                bgs += bgargv.strip().split(",")
            else:
                bgs.append(bgargv)

        for bgfile in bgs:
            if bgfile.endswith(".csv"):
                df = pd.read_csv(bgfile, dtype={"chrom": str})
                if "pos" in df.columns:
                    df["chrpos"] = df["chrom"] + "_" + df["pos"].astype(str)
                else:  # compatible with previous version
                    df["chrpos"] = df["chrom"] + "_" + df["pos2"].astype(str)
                df1 = df[["chrpos", "convs"]]
                df1.set_index("chrpos", inplace=True)
                for key1 in df1.index.to_list():
                    outdict[key1] = 1
            elif bgfile.endswith(".vcf"):
                bcf_in = pysam.VariantFile(bgfile)
                for rec in bcf_in.fetch():
                    try:
                        chrom, pos = rec.chrom, rec.pos
                        chr_pos = chrom + "_" + str(pos - 1)
                        outdict[chr_pos] = 1
                    except (ValueError, KeyError):
                        continue
                bcf_in.close()
            else:
                raise ValueError(
                    "Background snp file format cannot be recognized! Only csv or vcf format."
                )
        return outdict

    @utils.add_log
    def write_sparse_matrix(self, df, matrix_dir):
        count_matrix = CountMatrix.from_dataframe(
            df, self.features, barcodes=self.cell_list
        )
        count_matrix.to_matrix_dir(matrix_dir)

    @utils.add_log
    def write_h5ad(self):
        self.adata.layers["total"] = self.adata.X.copy()
        a1 = sc.read_10x_mtx(
            self.dir_labeled,
            var_names="gene_ids",
        )
        a1 = a1[self.adata.obs.index, self.adata.var.index]
        self.adata.layers["labeled"] = a1.X.copy()
        a2 = sc.read_10x_mtx(
            self.dir_unlabeled,
            var_names="gene_ids",
        )
        a2 = a2[self.adata.obs.index, self.adata.var.index]
        self.adata.layers["unlabeled"] = a2.X.copy()

        self.tor_stat()
        self.adata.write(self.h5ad)

    @utils.add_log
    def tor_stat(self):
        gene_count = self.adata.layers["total"].sum(axis=0)
        self.adata.var["total_counts"] = pd.DataFrame(gene_count).loc[0].to_numpy()
        self.adata = self.adata[:, self.adata.var["total_counts"] > 0]

        cell_ntr = self.adata.layers["labeled"].sum(axis=1) / self.adata.layers[
            "total"
        ].sum(axis=1)
        gene_ntr = self.adata.layers["labeled"].sum(axis=0) / self.adata.layers[
            "total"
        ].sum(axis=0)
        self.adata.obs["TOR"] = cell_ntr
        self.adata.var["TOR"] = gene_ntr.T

    @utils.add_log
    def tor_plot(self):
        tsne = pd.read_csv(self.tsne, sep="\t", index_col=0)
        tsne = tsne.loc[self.adata.obs.index]
        self.adata.obs["tSNE_1"] = tsne["tSNE_1"]
        self.adata.obs["tSNE_2"] = tsne["tSNE_2"]

        tsne_tor = Tsne_plot(
            self.adata.obs.sort_values(by="TOR"), "TOR", discrete=False
        )
        tsne_tor.set_color_scale("PuRd")
        self.add_data(tsne_tor=tsne_tor.get_plotly_div())

        vln_gene = Violin_plot(
            self.adata.var[self.adata.var["total_counts"] >= 10]["TOR"],
            "gene",
            color="#1f77b4",
        ).get_plotly_div()
        self.add_data(violin_gene=vln_gene)
        vln_cell = Violin_plot(
            self.adata.obs["TOR"], "cell", color="#ff7f0e"
        ).get_plotly_div()
        self.add_data(violin_cell=vln_cell)

    @utils.add_log
    def add_help(self):
        self.add_help_content(
            name="TOR:",
            content="(RNA turn-over rate) Fraction of labeled transcripts per gene or cell.",
        )

    def run_cmd(self, cmd):
        subprocess.call(" ".join(cmd), shell=True)

    @utils.add_log
    def clean_tmp(self):
        cmd = f"rm -rf {self.tmp_dir}"
        self.debug_subprocess_call(cmd)


@utils.add_log
def replacement(args):
    if args.control:
        return
    with Replacement(args, display_title="Labeled") as runner:
        runner.run()


def get_opts_replacement(parser, sub_program):
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"])
    parser.add_argument(
        "--control",
        action="store_true",
        help="For control samples to generate backgroup snp files and skip replacement",
    )
    parser.add_argument(
        "--replacementMem",
        default=50,
        type=int,
        help="Set replacement memory.",
    )
    if sub_program:
        parser.add_argument(
            "--bam",
            help="BAM file from the conversion step (if the temp dir was removed)",
            required=False,
        )
        parser.add_argument(
            "--cell", help="barcode cell list(from filtered matrix dir)", required=True
        )
        parser.add_argument(
            "--gene", help="gene list(from filtered matrix dir)", required=True
        )
        parser.add_argument(
            "--bg",
            nargs="+",
            required=False,
            help="background snp file, csv or vcf format",
        )
        parser.add_argument(
            "--tsne", help="tsne file from analysis step", required=True
        )
        parser = s_common(parser)
    return parser
