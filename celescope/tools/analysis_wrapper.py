import scanpy as sc
import numpy as np
import pandas as pd

from celescope.tools import utils
from celescope.__init__ import HELP_DICT
from celescope.rna.mkref import Mkref_rna
from celescope.tools.step import Step, s_common

# markers adjust p_value
PVAL_CUTOFF = 0.05
# scanpy mitochonrial variable name
MITO_VAR = 'mito'
NORMALIZED_LAYER = 'normalised'
RESOLUTION = 1.2
N_PCS = 25
MITO_GENE_PERCENT_LIST = [5, 10, 15, 20, 50]
# output marker top n in html
MARKER_TOP_N = 50
# marker sort by
MARKER_SORT_BY = 'p_val_adj'


def read_tsne(tsne_file):
    df = pd.read_csv(tsne_file, sep='\t')
    # compatible with old version
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={'Unnamed: 0': 'barcode'}, inplace=True)
        df = df.set_index('barcode')
    return df

def format_df_marker(df_marker):

    avg_logfc_col = "avg_log2FC"  # seurat 4
    if "avg_logFC" in df_marker.columns:  # seurat 2.3.4
        avg_logfc_col = "avg_logFC"
    df_marker = df_marker.loc[:,
                                ["cluster", "gene", avg_logfc_col, "pct.1", "pct.2", "p_val_adj"]
                                ]
    df_marker["cluster"] = df_marker["cluster"].apply(lambda x: f"cluster {x}")
    df_marker = df_marker[df_marker["p_val_adj"] < PVAL_CUTOFF]

    return df_marker

def get_opts_analysis(parser, sub_program):
    
    parser.add_argument('--genomeDir', help=HELP_DICT['genomeDir'], required=True)
    if sub_program:
        parser.add_argument(
            '--matrix_file',
            help='Required. Matrix_10X directory from step count.',
            required=True,
        )
        parser = s_common(parser)


class Scanpy_wrapper(Step):

    def __init__(self, args, display_title=None):

        super().__init__(args, display_title=display_title)

        # data
        self.adata = sc.read_10x_mtx(
            args.matrix_file,  
            var_names='gene_symbols',
        )
        self.mt_gene_list = Mkref_rna.parse_genomeDir(args.genomeDir)['mt_gene_list']

        # out
        self.df_marker_file = f'{self.out_prefix}_markers.tsv'
        self.df_marker_raw_file = f'{self.out_prefix}_markers_raw.tsv'
        self.df_tsne_file = f'{self.out_prefix}_tsne_coord.tsv'
        self.h5ad_file = f'{self.out_prefix}.h5ad'

   
    @utils.add_log
    def calculate_qc_metrics(self):

        if self.mt_gene_list:
            mito_genes, _ = utils.read_one_col(self.mt_gene_list)
            self.adata.var[MITO_VAR] = self.adata.var_names.map(lambda x:True if x in mito_genes else False)
            # if not astype(bool), it will be type object and raise an error
            # https://github.com/theislab/anndata/issues/504
            self.adata.var[MITO_VAR] = self.adata.var[MITO_VAR].astype(bool)
        else:
            self.adata.var[MITO_VAR] = self.adata.var_names.str.upper().str.startswith('MT-')

        sc.pp.calculate_qc_metrics(
            self.adata,
            qc_vars=[MITO_VAR], 
            percent_top=None,
            use_raw=False,
            log1p=False, 
            inplace=True
        )

    @utils.add_log
    def write_mito_stats(self):

        mt_pct_var = f'pct_counts_{MITO_VAR}'
        total_cell_number = self.adata.n_obs

        for mito_gene_percent in MITO_GENE_PERCENT_LIST:
            cell_number = sum(self.adata.obs[mt_pct_var] > mito_gene_percent)
            fraction = round(cell_number / total_cell_number * 100, 2)
            self.add_metric(
                name=f'Fraction of cells have mito gene percent>{mito_gene_percent}%',
                value=f'{fraction}%',
            )


    @utils.add_log
    def normalize(self):
        """
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
            n_pcs=N_PCS,
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
            n_pcs=N_PCS,
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
    def find_marker_genes(self):
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
        '''
        write only p_val_adj < PVAL_CUTOFF to avoid too many markers
        '''
        df_markers = sc.get.rank_genes_groups_df(self.adata, group=None, pval_cutoff=PVAL_CUTOFF)
        df_markers = df_markers[df_markers['logfoldchanges'].notna()]
        markers_name_dict = {
            'group':'cluster',
            'names':'gene',
            'logfoldchanges':'avg_log2FC',
            'pvals':'p_val',
            'pvals_adj':'p_val_adj',
            'pct_nz_group':'pct.1',
            'pct_nz_reference':'pct.2'
         }
        df_markers = df_markers.rename(markers_name_dict,axis='columns')
        df_markers['cluster'] = df_markers['cluster'].map(lambda x : int(x)+1)
        df_markers = df_markers.loc[df_markers['p_val_adj'] < PVAL_CUTOFF, ]
        df_markers.to_csv(self.df_marker_raw_file, index=None, sep='\t')

        df_markers_filter = df_markers.loc[df_markers['avg_log2FC'] > 0].sort_values(MARKER_SORT_BY, ascending=False).groupby('cluster').head(50)
        df_markers_filter = df_markers_filter.round({
            'avg_log2FC':3,
            'pct.1':3,
            'pct.2':3,
        })
        df_markers_filter.to_csv(self.df_marker_file, index=None, sep='\t')


    @utils.add_log
    def write_tsne(self):
        df_tsne = self.adata.obsm.to_df()[['X_tsne1','X_tsne2']]
        df_tsne['cluster']=self.adata.obs.cluster
        df_tsne['Gene_Counts']=self.adata.obs.n_genes_by_counts
        tsne_name_dict={'X_tsne1':'tSNE_1', 'X_tsne2':'tSNE_2'}
        df_tsne = df_tsne.rename(tsne_name_dict,axis='columns')
        df_tsne['cluster'] = df_tsne['cluster'].map(lambda x:int(x)+1)
        df_tsne.to_csv(self.df_tsne_file, sep='\t')

    @utils.add_log
    def write_h5ad(self):
        self.adata.write(self.h5ad_file)
    

    @utils.add_log
    def run(self):

        self.calculate_qc_metrics()
        self.write_mito_stats()
        self.normalize()
        self.hvg()
        self.scale()
        self.pca()
        self.neighbors()
        self.tsne()
        self.umap()
        self.leiden()
        self.find_marker_genes()

        self.write_markers()
        self.write_tsne()
        self.write_h5ad()


    def get_df(self):
        """
        return df_tsne, df_marker
        """
        df_tsne = read_tsne(self.df_tsne_file)
        df_marker = pd.read_csv(self.df_marker_file, sep="\t")
        df_marker = format_df_marker(df_marker)
        return df_tsne, df_marker


def get_opts_analysis_match(parser, sub_program):
    """
    Do not perform analysis. Only read data from scRNA-seq match_dir.
    """
    if sub_program:
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'])
        parser.add_argument("--tsne_file", help=HELP_DICT['tsne_file'])
        parser.add_argument("--df_marker_file", help=HELP_DICT['df_marker_file'])

        parser = s_common(parser)

class Report_runner(Step):

    def __init__(self, args, display_title=None):
    
        super().__init__(args, display_title=display_title)

    def add_marker_help(self):
        self.add_help_content(
            name='Marker Genes by Cluster',
            content='differential expression analysis based on the non-parameteric Wilcoxon rank sum test'
        )
        self.add_help_content(
            name='avg_log2FC',
            content='log fold-change of the average expression between the cluster and the rest of the sample'
        )
        self.add_help_content(
            name='pct.1',
            content='The percentage of cells where the gene is detected in the cluster'
        )
        self.add_help_content(
            name='pct.2',
            content='The percentage of cells where the gene is detected in the rest of the sample'
        )
        self.add_help_content(
            name='p_val_adj',
            content='Adjusted p-value, based on bonferroni correction using all genes in the dataset'
        )

    @staticmethod
    def get_df_file(match_dir):
        """
        return df_tsne_file, df_marker_file
        """
        match_dict = utils.parse_match_dir(match_dir)
        df_tsne_file = match_dict['tsne_coord']
        df_marker_file = match_dict.get('markers', None)
        return df_tsne_file, df_marker_file

    def get_df(self):
        """
        return df_tsne, df_marker
        """
        if utils.check_arg_not_none(self.args, 'match_dir'):
            df_tsne_file, df_marker_file = self.get_df_file(self.args.match_dir)
        elif utils.check_arg_not_none(self.args, 'tsne_file'):
            df_tsne_file = self.args.tsne_file
            df_marker_file = self.args.df_marker_file
        else:
            raise ValueError('match_dir or tsne_file must be specified')
        df_tsne = read_tsne(df_tsne_file)
        if df_marker_file:
            df_marker = pd.read_csv(df_marker_file, sep="\t")
            df_marker = format_df_marker(df_marker)
        else:
            df_marker = None
        return df_tsne, df_marker

    def run(self):
        pass
