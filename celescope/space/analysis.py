import sys

import numpy as np
import pandas as pd
from celescope.tools.step import s_common
import scanpy as sc
import scanpy.experimental as skp
import matplotlib.pyplot as plt
from matplotlib import colors
from celescope.tools.utils import add_log
from celescope.tools.plotly_plot import StaticPlot
from celescope.tools.analysis_wrapper import (
    Analysis as Tools_analysis,
    format_df_marker,
    NORMALIZED_LAYER,
    COUNTS_LAYER,
)
from celescope.__init__ import HELP_DICT
from celescope.space.utils import Spatial, convert_10x_h5


def hires_nocrop_spatial(adata, **kwargs):
    """
    sc.pl.spatial wrapper: show full hires image (no crop)
    """
    spatial = adata.uns["spatial"]["library0"]
    img = spatial["images"]["hires"]
    h, w = img.shape[:2]

    # scale factor (align coords)
    scale = spatial["scalefactors"].get("tissue_hires_scalef", 1.0)

    return sc.pl.spatial(
        adata,
        img_key="hires",
        show=False,
        crop_coord=(0, w / scale, 0, h / scale),
        **kwargs,
    )


class Analysis(Tools_analysis):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.filtered_counts_png = f"{self.outdir}/counts.png"
        self.raw_counts_png = f"{self.outdir}/raw_counts.png"
        self.cluster_png = f"{self.outdir}/cluster.png"
        self.raw_h5 = f"{self.outdir}/raw_feature_bc_matrix.h5"
        self.filtered_h5 = f"{self.outdir}/filtered_feature_bc_matrix.h5"
        self.spatial_dir = f"{self.outdir}/spatial"
        self.bad_channel_row_col_png = f"{self.outdir}/bad_channel_row_col.png"
        self.bad_channel_spatial_png = f"{self.outdir}/bad_channel_spatial.png"
        self.outs = [
            self.df_marker_raw_file,
            self.h5ad_file,
            self.filtered_counts_png,
            self.cluster_png,
            self.raw_h5,
            self.filtered_h5,
            self.spatial_dir,
        ]

    @add_log
    def read_data(self):
        self.adata = sc.read_visium(self.outdir, count_file="raw_feature_bc_matrix.h5")
        self.adata.layers[COUNTS_LAYER] = self.adata.X.copy()
        sc.pp.filter_genes(self.adata, min_cells=3)
        sc.pp.filter_cells(self.adata, min_genes=self.args.min_genes)
        self.adata.layers[NORMALIZED_LAYER] = self.adata.X.copy()

    @add_log
    def convert_h5(self):
        convert_10x_h5(self.args.raw, self.raw_h5)
        convert_10x_h5(self.args.filtered, self.filtered_h5)

    @add_log
    def get_filtered_data(self):
        adata_filtered = self.adata[self.adata.obs["in_tissue"] == 1, :].copy()
        self.adata = adata_filtered

    @add_log
    def run(self):
        self.convert_h5()
        spatial = Spatial(self.args.spatial)
        spatial.output_spatial(self.spatial_dir)
        self.read_data()
        self.calculate_qc_metrics()
        self.add_count_plot(self.raw_counts_png)
        self.get_filtered_data()
        if self.args.debug:
            self.bad_channel_analysis()
        self.add_count_plot(self.filtered_counts_png)
        self.add_data(plotly_count=StaticPlot(self.filtered_counts_png).get_div())
        self.write_mito_stats()
        # similar to Seurat sctransform
        skp.pp.normalize_pearson_residuals(self.adata)
        skp.pp.highly_variable_genes(
            self.adata, flavor="pearson_residuals", n_top_genes=2000
        )
        if np.isnan(self.adata.X).any():
            self.adata.X[np.isnan(self.adata.X)] = 0
        self.pca()
        self.neighbors()
        self.umap()
        self.leiden(resolution=0.3)
        self.add_cluster_plot()

        sc.pp.normalize_total(
            self.adata, target_sum=1e4, inplace=True, layer=NORMALIZED_LAYER
        )
        sc.pp.log1p(self.adata, layer=NORMALIZED_LAYER)
        n_clusters = len(self.adata.obs["cluster"].unique())
        if n_clusters > 1:
            sc.tl.rank_genes_groups(
                self.adata,
                "cluster",
                reference="rest",
                pts=True,
                use_raw=False,
                layer=NORMALIZED_LAYER,
                method="wilcoxon",
            )
            self.write_markers()
            self.add_marker_to_html()
        else:
            sys.stderr.write(
                "Warning: Only one cluster found. Skipping rank_genes_groups."
            )
        self.write_h5ad()
        spatial.rename_tissue_positions_csv(self.spatial_dir)

    @add_log
    def bad_channel_analysis(self):
        """
        Detect bad channels (rows/columns) with abnormally low UMI/gene counts
        and output images and metrics when args.debug is True.
        """
        X_dense = self.adata.layers[COUNTS_LAYER]
        self.adata.obs["total_UMI"] = np.array(X_dense.sum(axis=1)).flatten()
        self.adata.obs["n_genes"] = np.array((X_dense > 0).sum(axis=1)).flatten()

        # Step 1: aggregate by array_row
        row_stats = []
        for row, group in self.adata.obs.groupby("array_row"):
            if len(group) < 10:
                continue
            row_stats.append(
                [
                    row,
                    group["total_UMI"].median(),
                    group["n_genes"].median(),
                ]
            )
        row_stats = pd.DataFrame(
            row_stats,
            columns=["array_row", "median_total_UMI", "median_n_genes"],
        )
        umi_thresh_row = row_stats["median_total_UMI"].median() * 0.3
        gene_thresh_row = row_stats["median_n_genes"].median() * 0.3
        row_stats["bad_channel"] = (row_stats["median_total_UMI"] < umi_thresh_row) | (
            row_stats["median_n_genes"] < gene_thresh_row
        )

        # Step 2: aggregate by array_col
        col_stats = []
        for col, group in self.adata.obs.groupby("array_col"):
            if len(group) < 10:
                continue
            col_stats.append(
                [
                    col,
                    group["total_UMI"].median(),
                    group["n_genes"].median(),
                ]
            )
        col_stats = pd.DataFrame(
            col_stats,
            columns=["array_col", "median_total_UMI", "median_n_genes"],
        )
        umi_thresh_col = col_stats["median_total_UMI"].median() * 0.3
        gene_thresh_col = col_stats["median_n_genes"].median() * 0.3
        col_stats["bad_channel"] = (col_stats["median_total_UMI"] < umi_thresh_col) | (
            col_stats["median_n_genes"] < gene_thresh_col
        )

        # Step 3: mark spots in bad channels
        self.adata.obs["bad_row_channel"] = self.adata.obs["array_row"].isin(
            row_stats.loc[row_stats["bad_channel"], "array_row"]
        )
        self.adata.obs["bad_col_channel"] = self.adata.obs["array_col"].isin(
            col_stats.loc[col_stats["bad_channel"], "array_col"]
        )
        self.adata.obs["bad_channel"] = (
            self.adata.obs["bad_row_channel"] | self.adata.obs["bad_col_channel"]
        )

        # Step 4: visualize
        plt.figure(figsize=(12, 4))
        plt.subplot(1, 2, 1)
        plt.bar(
            row_stats["array_row"],
            row_stats["median_total_UMI"],
            color=np.where(row_stats["bad_channel"], "red", "gray"),
        )
        plt.xlabel("array_row")
        plt.ylabel("Median total UMI")
        plt.title("Row (X) UMI")

        plt.subplot(1, 2, 2)
        plt.bar(
            col_stats["array_col"],
            col_stats["median_total_UMI"],
            color=np.where(col_stats["bad_channel"], "red", "gray"),
        )
        plt.xlabel("array_col")
        plt.ylabel("Median total UMI")
        plt.title("Column (Y) UMI")
        plt.tight_layout()
        plt.savefig(self.bad_channel_row_col_png, dpi=300, bbox_inches="tight")
        plt.close()

        x = self.adata.obs["array_col"]
        y = self.adata.obs["array_row"]
        spot_colors = self.adata.obs["bad_channel"].map({True: "red", False: "gray"})
        plt.figure(figsize=(8, 8))
        plt.scatter(x, y, c=spot_colors, s=10, edgecolor="k", linewidth=0.2)
        plt.gca().invert_yaxis()
        plt.xlabel("array_col")
        plt.ylabel("array_row")
        plt.title("Spatial plot of bad channels")
        plt.axis("equal")
        plt.savefig(self.bad_channel_spatial_png, dpi=300, bbox_inches="tight")
        plt.close()

        # Step 5: output stats
        bad_channel_percent = (
            self.adata.obs["bad_channel"].sum() / self.adata.n_obs * 100
        )
        bad_row_str = f"{row_stats['bad_channel'].sum()} / {len(row_stats)}"
        bad_col_str = f"{col_stats['bad_channel'].sum()} / {len(col_stats)}"
        bad_spot_str = (
            f"{self.adata.obs['bad_channel'].sum()} / {self.adata.n_obs} "
            f"({bad_channel_percent:.2f}%)"
        )
        self.add_metric(
            name="bad row channels",
            value=bad_row_str,
            help_info="Number of row channels flagged as bad vs total rows",
            show=False,
        )
        self.add_metric(
            name="bad column channels",
            value=bad_col_str,
            help_info="Number of column channels flagged as bad vs total columns",
            show=False,
        )
        self.add_metric(
            name="bad channel spots",
            value=bad_spot_str,
            help_info="Number of spots flagged as bad channels vs total spots",
            show=False,
        )

        row_stats.to_csv(f"{self.outdir}/bad_channel_row_stats.csv", index=False)
        col_stats.to_csv(f"{self.outdir}/bad_channel_col_stats.csv", index=False)

    @add_log
    def add_count_plot(self, plot_path):
        plt.figure(figsize=(8, 8))
        hires_nocrop_spatial(
            self.adata,
            color=["total_counts"],
            color_map="Reds",
            size=1.5,
            alpha=0.5,
            norm=colors.LogNorm(vmin=1),
        )
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

    @add_log
    def add_cluster_plot(self):
        plt.figure(figsize=(8, 8))
        hires_nocrop_spatial(self.adata, color=["cluster"], size=1.5)
        plt.savefig(self.cluster_png, dpi=300, bbox_inches="tight")
        plt.close()
        self.add_data(plotly_cluster=StaticPlot(self.cluster_png).get_div())

    @add_log
    def add_marker_to_html(self):
        df_marker = pd.read_csv(self.df_marker_file, sep="\t")
        df_marker = format_df_marker(df_marker)
        script = """
<script>

    $(document).ready(function () {
            var table = $('#marker_genes').DataTable({
                dom: 'Bfrtip',
                buttons: ['excel']
            });
            var indexOfMyCol = 0 ;
            var collator = new Intl.Collator(undefined, {numeric: true, sensitivity: 'base'});
    $("#marker_genes thead th").each( function ( i ) {
        if (i==indexOfMyCol){

          var select = $('<select><option value=""></option></select>')
            .appendTo( $(this).empty() )
            .on( 'change', function () {
                var pattern = ""
                if ($(this).val()!="") {
                    pattern= pattern="^"+$(this).val() +"$"
                }
                table.column( i )
                .search(input=pattern, regex=true, smart=false)
                .draw();
            } );
 
        table.column( i).data().unique().sort(collator.compare).each( function ( d, j ) {
            select.append( '<option value="'+d+'">'+d+'</option>' )
        } );
    }
    } );
    });
</script>
"""
        self.add_table(
            title="Marker Genes by Cluster",
            table_id="marker_genes",
            df=df_marker,
            script=script,
        )


def analysis(args):
    with Analysis(args) as runner:
        runner.run()


def get_opts_analysis(parser, sub_program):
    parser.add_argument(
        "--min_genes",
        help="Minimum number of genes expressed required for a cell to pass filtering.",
        type=int,
        default=1,
    )
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"], required=True)
    if sub_program:
        parser.add_argument("--raw", help="raw matrix", required=True)
        parser.add_argument("--filtered", help="filtered matrix", required=True)
        parser.add_argument("--spatial", help="spatial directory", required=True)
        parser = s_common(parser)
