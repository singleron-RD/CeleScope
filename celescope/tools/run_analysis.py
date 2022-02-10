import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from celescope.tools import utils


class Scanpy():
    def __init__(self,matrix_file,outdir,sample,mt_gene_list=None,save_h5ad=False) -> None:
        
        DIMS = 20
        RESOLUTION = 0.6
        N_FEATURES = 20000
        
        self.matrix_file = matrix_file
        self.outdir = outdir
        self.sample = sample
        self.mt_gene_list = mt_gene_list
        self.save_h5ad = save_h5ad
    
    def filter(self,adata, min_genes=200, max_genes=5000, pct_counts_mito=50, *args, **kwargs):
        """
        Wrapper function for sc.pp.filter_cells() and sc.pp.filter_genes(), mainly
        for supporting arbitrary filtering
        """
        layer = 'counts' if 'counts' in adata.layers.keys() else None
        pct_top = []

        sc.pp.filter_cells(adata, min_genes=min_genes, inplace=True, copy=False)
        sc.pp.filter_cells(adata, max_genes=max_genes, inplace=True, copy=False)
        # sc.pp.filter_cells(adata, min_counts=0, inplace=True, copy=False)
        # sc.pp.filter_cells(adata, max_counts=30000, inplace=True, copy=False)
        sc.pp.filter_genes(adata, min_cells=5, inplace=True, copy=False)

        sc.pp.calculate_qc_metrics(
            adata,
            expr_type='counts',
            var_type='genes',
            qc_vars=['mt'], 
            percent_top=pct_top,
            layer=layer,
            use_raw=False
            log1p=False, 
            inplace=True
        )

        if 'mt' in adata.var:
            adata = adata[adata.obs['pct_counts_mt'] < pct_counts_mito, :]

        return adata

    def normalize(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.pp.normalize_per_cell() and sc.pp.log1p(), mainly
        for supporting different ways of saving raw data.
        """
        adata.layers['filtered'] = adata.X
        sc.pp.normalize_total(
            adata,
            target_sum=1e4,
            exclude_highly_expressed=False,
            max_fraction=0.05,
            key_added=None,
            layer=None,
            layer_norm=None,
            inplace=True,
            copy=False,
        )
        sc.pp.log1p(
            adata,
            base=None,
            copy=False,
            chunked=None,
            chunk_size=None,
            layer=None,
            obsm=None
        )

        return adata

    def hvg(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.highly_variable_genes()
        """
        sc.pp.highly_variable_genes(
            adata,
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

        return adata

    def scale(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.pp.scale
        """
        adata.layers['normalised'] = adata.X
        sc.pp.scale(
            adata,
            zero_center=True,
            max_value=10,
            copy=False,
            layer=None,
            obsm=None
        )

        return adata

    def pca(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.pp.pca
        """
        sc.pp.pca(
            adata,
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

        return adata

    def harmony(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.external.pp.harmony_integrate
        """
        sc.external.pp.harmony_integrate(
            adata,
            key='Sample ID',
            basis='X_pca',
            adjusted_basis='X_pca',
        )

        return adata

    def neighbors(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.pp.neighbors(), for supporting multiple n_neighbors
        """
        sc.pp.neighbors(
            adata,
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

        return adata

    def tsne(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.tl.tsne, for supporting named slot of tsne embeddings
        """
        sc.tl.tsne(
            adata,
            n_pcs=50,
            use_rep=None,
            perplexity=30,
            early_exaggeration=12,
            learning_rate=1000,
            random_state=0,
            use_fast_tsne=False,
            n_jobs=None,
            copy=False,
            metric='euclidean'
        )

        return adata

    def umap(self,adata, *args, **kwargs):
        """
        Wrapper function for sc.tl.umap, for supporting named slot of umap embeddings
        """
        sc.tl.umap(
            adata,
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

        return adata

    def leiden(self,adata, resolution=1, *args, **kwargs):
        """
        Wrapper function for sc.tl.leiden, for supporting multiple resolutions.
        """
        sc.tl.leiden(
            adata,
            resolution=resolution,
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

        return adata

    def run(self):
        adata = sc.read_10x_mtx(
                    self.matrix_file,  
                    var_names='gene_symbols',
                    cache=True)
        
        all_genes = adata.var['gene_ids'].index.to_list()

        if self.mt_gene_list != None:
            mt,_ = read_one_col(self.mt_gene_list)
            mito_genes = list(set(mt).intersection(set(all_genes)))
            adata.var['mt'] =adata.var_names.map(lambda x:True if x in mito_genes else False)
        else:
            adata.var['mt'] = adata.var_names.str.startswith('MT-')

        adata = self.filter(adata)

        #stat.txt
        total_cell = len(adata.obs)
        percent_list = [5,10,15,20,50]
        mito_dict = {key:sum(adata.obs['pct_counts_mt'] > key)/total_cell for key in percent_list}
        mito_df = pd.DataFrame.from_dict(mito_dict,orient='index',columns=['cell_percent'])
        mito_df = mito_df.reset_index().rename(columns={'index':'percent_mito'})
        mito_df['cell_percent'] = round(mito_df['cell_percent']*100,2).astype(str)
        mito_df['cell_percent'] += '%'
        mito_df['percent_mito'] = round(mito_df['percent_mito'],2).astype(str)
        mito_df['percent_mito'] += '%'
        mito_df['percent_mito'] = "Fraction of cells have mito gene percent>" + mito_df['percent_mito']
        mito_df.to_csv(f'{self.outdir}/stat.txt',sep=':',index=None,header=False)

        #analysis
        adata = self.normalize(adata)
        adata = self.hvg(adata)
        adata = self.scale(adata)
        adata = self.pca(adata)
        #if batch_remove:
        #    adata = self.harmony(adata)
        adata = self.neighbors(adata)
        adata = self.tsne(adata)
        adata = self.umap(adata)
        adata = self.leiden(adata)

        sc.tl.rank_genes_groups(adata,'leiden',use_highly_variable=True,method='wilcoxon',
                                corr_method='bonferroni',n_genes=None,pts=True)

        #report
        #markers.tsv
        df_markers = sc.get.rank_genes_groups_df(adata,group=None)
        df_markers = df_markers[df_markers['logfoldchanges'].notna()]
        markers_name_dict={'group':'cluster','names':'gene','logfoldchanges':'avg_log2FC',
                           'pvals':'p_val','pvalse_adj':'p_val_adj',
                           'pct_nz_group':'pct.1','pct_nz_reference':'pct.2'}
        df_markers = df_markers.rename(markers_name_dict,axis='columns')
        df_markers['cluster'] = df_markers['cluster'].map(lambda x : int(x)+1)
        df_markers.to_csv(f'{self.outdir}/{self.sample}_markers.tsv',index=None,sep='\t')

        #tsne_coord.tsv
        df_tsne = adata.obsm.to_df()[['X_tsne1','X_tsne2']]
        df_tsne['cluster']=adata.obs.louvain
        df_tsne['Gene_Counts']=adata.obs.n_genes_by_counts
        tsne_name_dict={'X_tsne1':'tSNE_1','X_tsne2':'tSNE_2'}
        df_tsne = df_tsne.rename(tsne_name_dict,axis='columns')
        df_tsne['cluster'] = df_tsne['cluster'].map(lambda x:int(x)+1)
        df_tsne.to_csv(f'{self.outdir}/{self.sample}_tsne_coord.tsv',sep='\t')

        #save h5ad file
        results_file = f'{self.outdir}/{self.sample}.h5ad'
        if self.save_h5ad != None:
            adata.write(results_file)