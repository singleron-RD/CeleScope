from celescope.tools.step import Step, s_common
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import colors
from skimage.color import gray2rgb
from celescope.tools.utils import add_log


class Analysis(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.counts_png = f"{self.outdir}/counts.png"
        self.cluster_png = f"{self.outdir}/cluster.png"
        self.outs = [self.counts_png, self.cluster_png]

    @add_log
    def run(self):
        adata = sc.read_visium(
            self.args.outs_dir, count_file="filtered_feature_bc_matrix.h5"
        )

        adata.var["mt"] = adata.var_names.str.startswith("mt-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

        lib_id = list(adata.uns["spatial"].keys())[0]
        img = adata.uns["spatial"][lib_id]["images"]["hires"]

        # 如果是灰度图
        if img.ndim == 2:
            adata.uns["spatial"][lib_id]["images"]["hires"] = gray2rgb(img)

        # 1. QC 图片
        sc.pl.spatial(
            adata,
            img_key="hires",
            color=["total_counts"],
            color_map="jet",
            size=1.5,
            norm=colors.LogNorm(vmin=1),
            show=False,
            save=None,
        )
        plt.savefig(self.counts_png, dpi=300, bbox_inches="tight")
        plt.close()

        sc.pp.filter_cells(adata, min_genes=self.args.min_genes)

        # 2. cluster 图片
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.8)

        plt.rcParams["figure.figsize"] = (8, 8)
        sc.pl.spatial(
            adata,
            color=["leiden"],
            img_key="hires",
            size=1.5,
            show=False,
            cmap="tab20",
            save=None,
        )
        plt.savefig(self.cluster_png, dpi=300, bbox_inches="tight")
        plt.close()


def analysis(args):
    with Analysis(args) as runner:
        runner.run()


def get_opts_analysis(parser, sub_program):
    parser.add_argument(
        "--min_genes",
        help="Minimum number of genes expressed required for a cell to pass filtering.",
        type=int,
        default=200,
    )
    if sub_program:
        parser.add_argument("--outs_dir", help="celescope outs dir", required=True)
        parser = s_common(parser)
