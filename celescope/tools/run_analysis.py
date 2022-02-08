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

        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

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
        sc.pp.normalize_total(adata,target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata,n_top_genes=N_FEATURES,min_mean=0.1,max_mean=8,min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)  #50?
        sc.tl.pca(adata, svd_solver='arpack',use_highly_variable=True)
        sc.pp.neighbors(adata,n_neighbors=15,n_pcs=DIMS)
        sc.tl.louvain(adata,resolution=RESOLUTION)
        sc.tl.tsne(adata,n_pcs=DIMS,learning_rate=1000) #learning_rate is 200 in seurat 
        sc.tl.umap(adata,min_dist=0.3)
        sc.tl.rank_genes_groups(adata,'louvain',use_highly_variable=True,method='wilcoxon',
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