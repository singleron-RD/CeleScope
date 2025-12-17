import pandas as pd
from celescope.tools.step import s_common
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import colors
from celescope.tools.utils import add_log
from celescope.tools.plotly_plot import StaticPlot
from celescope.tools.analysis_wrapper import Scanpy_wrapper, format_df_marker
from celescope.__init__ import HELP_DICT
from celescope.space.utils import Spatial, convert_10x_h5


class Analysis(Scanpy_wrapper):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.filtered_counts_png = f"{self.outdir}/counts.png"
        self.raw_counts_png = f"{self.outdir}/raw_counts.png"
        self.cluster_png = f"{self.outdir}/cluster.png"
        self.raw_h5 = f"{self.outdir}/raw_feature_bc_matrix.h5"
        self.filtered_h5 = f"{self.outdir}/filtered_feature_bc_matrix.h5"
        self.spatial_dir = f"{self.outdir}/spatial"
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

    @add_log
    def filter_cells(self):
        sc.pp.filter_cells(self.adata, min_genes=self.args.min_genes, inplace=True)

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
        self.add_count_plot(self.filtered_counts_png)
        self.add_data(plotly_count=StaticPlot(self.filtered_counts_png).get_div())
        self.write_mito_stats()
        self.filter_cells()
        self.normalize()
        self.hvg()
        self.scale()
        self.pca()
        self.neighbors()
        self.umap()
        self.leiden(resolution=0.8)
        self.add_cluster_plot()
        self.find_marker_genes()
        self.write_markers()
        self.write_h5ad()
        self.add_marker_to_html()
        spatial.rename_tissue_positions_csv(self.spatial_dir)

    @add_log
    def add_count_plot(self, plot_path):
        plt.figure(figsize=(8, 8))
        sc.pl.spatial(
            self.adata,
            img_key="hires",
            color=["total_counts"],
            color_map="Reds",
            size=1.5,
            norm=colors.LogNorm(vmin=1),
            show=False,
            save=None,
        )
        plt.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close()

    @add_log
    def add_cluster_plot(self):
        plt.figure(figsize=(8, 8))
        sc.pl.spatial(
            self.adata,
            color=["cluster"],
            img_key="hires",
            size=1.5,
            show=False,
            save=None,
        )
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
        default=200,
    )
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"], required=True)
    if sub_program:
        parser.add_argument("--raw", help="raw matrix", required=True)
        parser.add_argument("--filtered", help="filtered matrix", required=True)
        parser.add_argument("--spatial", help="spatial directory", required=True)
        parser = s_common(parser)
