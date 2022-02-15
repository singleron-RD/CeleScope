import subprocess
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd

from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.rna.mkref import Mkref_rna
from celescope.tools.step import Step, s_common
from celescope.tools.run_analysis import Scanpy


MITO_VAR = 'mito'
NORMALIZED_LAYER = 'normalised'
FILTERED_LAYER = 'filtered'

MIN_GENES = 200
MAX_GENES = 5000
MIN_CELLS = 5
RESOLUTION = 1.2



def get_opts_analysis(parser, sub_program):
    
    parser.add_argument('--genomeDir', help=HELP_DICT['genomeDir'], required=True)
    if sub_program:
        parser.add_argument(
            '--matrix_file',
            help='Required. Matrix_10X directory from step count.',
            required=True,
        )
        parser = s_common(parser)


def get_opts_analysis_match(parser, sub_program):
    if sub_program:
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--tsne_file", help="match_dir t-SNE coord file. Not required when `--match_dir` is provided.")
        parser.add_argument("--df_marker_file", help="match_dir df_marker_file. Not required when `--match_dir` is provided.")

        parser = s_common(parser)


class Scanpy_wrapper(Step):

    def __init__(self, args, display_title=None):

        super().__init__(args, display_title=display_title)

        # data
        self.adata = sc.read_10x_mtx(
            args.matrix_file,  
            var_names='gene_symbols',
            cache=True)
        self.mt_gene_list = Mkref_rna.parse_genomeDir(args.genomeDir)['mt_gene_list']

        # out
        self.marker_file = f'{self.out_prefix}_marker.tsv'
        self.tsne_file = f'{self.out_prefix}_tsne.tsv'
        self.h5ad_file = f'{self.out_prefix}.h5ad'


    @utils.add_log
    def filter(self):
        """
        sc.pp.filter_cells() and sc.pp.filter_genes()
        """

        sc.pp.filter_cells(self.adata, min_genes=MIN_GENES,  max_genes=MAX_GENES, inplace=True, copy=False)
        sc.pp.filter_genes(self.adata, min_cells=MIN_CELLS, inplace=True, copy=False)
        self.adata.layers[FILTERED_LAYER] = self.adata.X
    
    @utils.add_log
    def calculate_qc_metrics(self):

        if self.mt_gene_list:
            mito_genes, _ = utils.read_one_col(self.mt_gene_list)
            self.adata.var[MITO_VAR] = self.adata.var_names.map(lambda x:True if x in mito_genes else False)
        else:
            self.adata.var[MITO_VAR] = self.adata.var_names.str.upper().str.startswith('MT-')

        sc.pp.calculate_qc_metrics(
            self.adata,
            qc_vars=[MITO_VAR], 
            use_raw=False,
            log1p=False, 
            inplace=True
        )

    @utils.add_log
    def write_mito_stats(self):
 
        #stat.txt
        total_cell = len(self.adata.obs)
        percent_list = [5,10,15,20,50]
        mito_dict = {
            key: sum(self.adata.obs['pct_counts_mt'] > key) / total_cell             
            for key in percent_list
        }

        for percent in percent_list:
            self.add_metric(
                name=f'Fraction of cells have mito gene percent>{percent}%',
                value=round(mito_dict[percent], 2),
            )

    @utils.add_log
    def normalize(self):
        """
        add `filtered` layer
        sc.pp.normalize_per_cell() and sc.pp.log1p()
        """

        sc.pp.normalize_total(
            self.adata,
            target_sum=1e4,
            inplace=True,
        )
        sc.pp.log1p(
            self.adata,
        )
        self.adata.layers[NORMALIZED_LAYER] = self.adata.X

    @utils.add_log
    def hvg(self):
        """
        Wrapper function for sc.highly_variable_genes()
        """
        sc.pp.highly_variable_genes(
            self.adata,
            layer=None,
            n_top_genes=None,
            min_disp=0.5,
            max_disp=np.inf,
            min_mean=0.0125,
            max_mean=3,
            span=0.3,
            n_bins=20,
            flavor='seurat',
            subset=False,
            inplace=True,
            batch_key=None,
            check_values=True
        )        

    @utils.add_log
    def scale(self):
        """
        Wrapper function for sc.pp.scale
        """
        sc.pp.scale(
            self.adata,
            zero_center=True,
            max_value=10,
            copy=False,
            layer=None,
            obsm=None
        )
        
    @utils.add_log
    def pca(self):
        """
        Wrapper function for sc.pp.pca
        """
        sc.pp.pca(
            self.adata,
            n_comps=50,
            zero_center=True,
            svd_solver='auto',
            random_state=0,
            return_info=False,
            use_highly_variable=True,
            dtype='float32',
            copy=False,
            chunked=False,
            chunk_size=None
        )       

    @utils.add_log
    def neighbors(self, ):
        """
        Wrapper function for sc.pp.neighbors(), for supporting multiple n_neighbors
        """
        sc.pp.neighbors(
            self.adata,
            n_neighbors=15,
            n_pcs=25,
            use_rep=None,
            knn=True,
            random_state=0,
            method='umap',
            metric='euclidean',
            key_added=None,
            copy=False
        )

    @utils.add_log
    def tsne(self):
        """
        Wrapper function for sc.tl.tsne, for supporting named slot of tsne embeddings
        """
        sc.tl.tsne(
            self.adata,
            n_pcs=50,
            use_rep=None,
            perplexity=30,
            early_exaggeration=12,
            learning_rate=1000,
            random_state=0,
            n_jobs=self.thread,
            copy=False,
        )

    @utils.add_log
    def umap(self, ):
        """
        Wrapper function for sc.tl.umap, for supporting named slot of umap embeddings
        """
        sc.tl.umap(
            self.adata,
            min_dist=0.5,
            spread=1.0,
            n_components=2,
            maxiter=None,
            alpha=1.0,
            gamma=1.0,
            negative_sample_rate=5,
            init_pos='spectral',
            random_state=0,
            a=None,
            b=None,
            copy=False,
            method='umap',
            neighbors_key=None
        )

        
    @utils.add_log
    def leiden(self):
        """
        Wrapper function for sc.tl.leiden
        """
        sc.tl.leiden(
            self.adata,
            resolution=RESOLUTION,
            restrict_to=None,
            random_state=0,
            key_added='cluster',
            adjacency=None,
            directed=True,
            use_weights=True,
            n_iterations=-1,
            partition_type=None,
            neighbors_key=None,
            obsp=None,
            copy=False
        )

    @utils.add_log
    def cluster(self):
        """
        Wrapper function for sc.tl.rank_genes_groups
        """
        sc.tl.rank_genes_groups(
            self.adata,
            "cluster",
            reference='rest',
            pts=True,
            method="wilcoxon",
            use_raw=False,
            layer=NORMALIZED_LAYER
        )

    @utils.add_log
    def write_markers(self):
        df_markers = sc.get.rank_genes_groups_df(self.adata,group=None)
        df_markers = df_markers[df_markers['logfoldchanges'].notna()]
        markers_name_dict={'group':'cluster','names':'gene','logfoldchanges':'avg_log2FC',
                           'pvals':'p_val','pvalse_adj':'p_val_adj',
                           'pct_nz_group':'pct.1','pct_nz_reference':'pct.2'}
        df_markers = df_markers.rename(markers_name_dict,axis='columns')
        df_markers['cluster'] = df_markers['cluster'].map(lambda x : int(x)+1)
        df_markers.to_csv(self.marker_file, index=None, sep='\t')

    @utils.add_log
    def write_tsne(self):
        df_tsne = self.adata.obsm.to_df()[['X_tsne1','X_tsne2']]
        df_tsne['cluster']=self.adata.obs.cluster
        df_tsne['Gene_Counts']=self.adata.obs.n_genes_by_counts
        tsne_name_dict={'X_tsne1':'tSNE_1','X_tsne2':'tSNE_2'}
        df_tsne = df_tsne.rename(tsne_name_dict,axis='columns')
        df_tsne['cluster'] = df_tsne['cluster'].map(lambda x:int(x)+1)
        df_tsne.to_csv(self.tsne_file, sep='\t')

    @utils.add_log
    def write_h5ad(self):
        self.adata.write(self.h5ad_file)
    

    @utils.add_log
    def run(self):

        self.filter()
        self.calculate_qc_metrics()
        self.add_mito_metric()
        self.normalize()
        self.hvg()
        self.scale()
        self.pca()
        self.neighbors()
        self.tsne()
        self.umap()
        self.leiden()
        self.cluster()

        self.write_markers()
        self.write_tsne()
        self.write_h5ad()

    def read_match_dir(self):
        """
        if match_dir is not self, should read match_dir at init
        if it is self, read at run_analysis - need to run seurat first
        """
        if self.match_dir:
            match_dict = utils.parse_match_dir(self.match_dir)
            tsne_df_file = match_dict['tsne_coord']
            self.df_tsne = pd.read_csv(tsne_df_file, sep="\t")
            self.df_tsne.rename(columns={"Unnamed: 0": "barcode"}, inplace=True)
            self.df_marker_file = match_dict['markers']
            self.read_format_df_marker()

