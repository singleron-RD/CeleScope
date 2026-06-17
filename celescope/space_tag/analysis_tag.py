import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import colors

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.tools.plotly_plot import StaticPlot
from celescope.__init__ import HELP_DICT


class Analysis_tag(Step):
    """
    ## Features
    - Combine spatial RNA-Seq clustering information with tag assignment.

    ## Output

    - `{sample}_spatial_tag.png` - Spatial plot showing tag abundance
    - `{sample}_cluster_spatial.png` - Spatial plot showing clusters

    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)
        self.umi_tag_file = args.umi_tag_file

        # out files
        self.spatial_tag_png = f"{self.outdir}/{self.sample}_spatial_tag.png"
        self.cluster_spatial_png = f"{self.outdir}/{self.sample}_cluster_spatial.png"

    def hires_nocrop_spatial(self, adata, **kwargs):
        """
        sc.pl.spatial wrapper: show full hires image (no crop)
        """
        spatial = adata.uns["spatial"]["library0"]
        img = spatial["images"]["hires"]
        h, w = img.shape[:2]
        scale = spatial["scalefactors"].get("tissue_hires_scalef", 1.0)

        return sc.pl.spatial(
            adata,
            img_key="hires",
            show=False,
            crop_coord=(0, w / scale, 0, h / scale),
            **kwargs,
        )

    def run(self):
        # Read rna.h5ad from match_dir/outs
        if not hasattr(self.args, "match_dir") or not self.args.match_dir:
            raise ValueError("--match_dir is required")

        rna_h5ad_path = os.path.join(self.args.match_dir, "outs/rna.h5ad")

        if not os.path.exists(rna_h5ad_path):
            raise FileNotFoundError(f"Could not find rna.h5ad: {rna_h5ad_path}")

        adata = sc.read(rna_h5ad_path)

        # Read umi_tag file (matrix format: index=barcode, columns=tags, values=UMI counts)
        df_umi_tag = pd.read_csv(self.umi_tag_file, sep="\t", index_col=0)

        # Get tag columns (all columns are tag names)
        tag_columns = df_umi_tag.columns.tolist()

        # Calculate total tag UMI per cell
        df_umi_tag["total_tag_umi"] = df_umi_tag.sum(axis=1)

        # Add tag columns to adata.obs
        for tag in tag_columns:
            adata.obs[tag] = df_umi_tag[tag].reindex(adata.obs_names).fillna(0)

        # Add total tag UMI column
        adata.obs["total_tag_umi"] = (
            df_umi_tag["total_tag_umi"].reindex(adata.obs_names).fillna(0)
        )

        # Save modified h5ad
        adata.write(f"{self.outdir}/{self.sample}_tag.h5ad")

        # Generate cluster spatial plot if cluster column exists
        if "cluster" in adata.obs.columns:
            plt.figure(figsize=(8, 8))
            self.hires_nocrop_spatial(adata, color=["cluster"], size=1.5)
            plt.savefig(self.cluster_spatial_png, dpi=300, bbox_inches="tight")
            plt.close()
            self.add_data(plotly_cluster=StaticPlot(self.cluster_spatial_png).get_div())

        # Generate spatial plot for total tag UMI
        plt.figure(figsize=(8, 8))
        self.hires_nocrop_spatial(
            adata,
            color=["total_tag_umi"],
            size=1.5,
            color_map="Reds",
            norm=colors.LogNorm(vmin=1),
        )
        plt.savefig(self.spatial_tag_png, dpi=300, bbox_inches="tight")
        plt.close()
        self.add_data(plotly_tag=StaticPlot(self.spatial_tag_png).get_div())


def get_opts_analysis_tag(parser, sub_program):
    if sub_program:
        parser.add_argument(
            "--umi_tag_file",
            help="`{sample}_umi_tag.tsv` from count_tag. ",
            required=True,
        )
        parser.add_argument("--match_dir", help=HELP_DICT["match_dir"], required=True)
        parser = s_common(parser)


@utils.add_log
def analysis_tag(args):
    with Analysis_tag(args, display_title="Analysis") as runner:
        runner.run()
