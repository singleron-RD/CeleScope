import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import glob
import pathlib
from scipy.sparse import coo_matrix
from scipy.io import mmwrite

from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT
from celescope.tools.plotly_plot import Tsne_dropdown_plot,Tsne_plot



def fun_filter(df):
    df['filter_num'] = np.mean(df,axis=1)+np.std(df,axis=1)
    for index,_ in df.iterrows():
        df.loc[index] = df.loc[index].mask(df.loc[index]<df.loc[index]['filter_num'],0)
        #df.loc[index].apply(lambda x: 0 if x>df.loc[index]['filter_num'] else x)
    df.drop(columns='filter_num',inplace=True)
    return df


def matrix_to_10X(df,outdir):
    """
    creat 10X format
    """
    pathlib.Path(outdir).mkdir(parents=True,exist_ok=True)
    matrix = df
    obs = pd.DataFrame(index=df.columns)
    var = pd.DataFrame(df.index,index=df.index,columns=['ADT'])
    mtx = coo_matrix(matrix, matrix.shape, int)
    mmwrite(f"{outdir}/matrix", mtx)
    obs.to_csv(f"{outdir}/barcodes.tsv", sep='\t', header=False)
    var.to_csv(f"{outdir}/genes.tsv", sep='\t', header=False)


def get_opts_analysis_cite(parser, sub_program):
    if sub_program:
        parser.add_argument('--match_dir', help=HELP_DICT['match_dir'], required=True)
        parser.add_argument('--citeseq_mtx', help='citeseq matrix .gz file', required=True)
        s_common(parser)


def analysis_cite(args):
    with Analysis_cite(args, display_title='Analysis') as runner:
        runner.run()


class Analysis_cite(Step):
    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        # input_file
        self.citeseq_mtx = args.citeseq_mtx
        self.tsne_coord = glob.glob(f'{args.match_dir}/*analysis*/*tsne_coord.tsv')[0]
        self.filter_feature_bc_matrix = glob.glob(f'{args.match_dir}/*count*/*filtered_feature_bc_matrix')[0]

        #out
        self.raw_tag_bc_matrix  = f'{args.outdir}/{args.sample}_raw_tag_bc_matrix'
        self.filtered_tag_bc_matrix  = f'{args.outdir}/{args.sample}_filtered_tag_bc_matrix'
        #self.filter_citeseq_mtx = f'{args.outdir}/{args.sample}_filter_citeseq.mtx.gz'
        


    def run(self):
        """
        
        """
        #filter raw
        df_raw = pd.read_csv(self.citeseq_mtx,sep="\t",index_col=0)
        df_filter = fun_filter(df_raw)
        #df_filter = df_raw.loc[:,~(df_raw==0).all(axis=0)]

        #merge feature and tag
        mdata_feature = sc.read_10x_mtx(self.filter_feature_bc_matrix)
        df_feature = mdata_feature.to_df().T
        df_merge_raw = pd.concat([df_feature,df_raw],axis=0)
        df_merge_filter = pd.concat([df_feature,df_filter],axis=0)

        #creat 10X
        matrix_to_10X(df_merge_raw,self.raw_tag_bc_matrix)
        matrix_to_10X(df_merge_filter,self.filtered_tag_bc_matrix)

        #normalize citeseq date 
        obs = pd.DataFrame(index=df_filter.columns)
        var = pd.DataFrame(df_filter.index,index=df_filter.index,columns=['ADT'])
        mdata_citeseq = ad.AnnData(np.array(df_filter.T),obs=obs,var=var)
        sc.pp.normalize_total(
                mdata_citeseq,
                target_sum=1e4,
                inplace=True,
            )
        sc.pp.log1p(mdata_citeseq)

        #merge tsne_coord and filter_citeseq date 
        df_tsne = pd.read_csv(self.tsne_coord,sep="\t")
        if 'Unnamed: 0' in df_tsne.columns:
            df_tsne.rename(columns={'Unnamed: 0': 'barcode'}, inplace=True)
            df_tsne = df_tsne.set_index('barcode')

        df_citeseq = mdata_citeseq.to_df()
        df_tsne = pd.concat([df_tsne,df_citeseq],axis=1)
        df_tsne.fillna(0,inplace=True)

        feature_name_list = df_citeseq.columns.to_list()

        #plot
        tsne_cluster = Tsne_plot(df_tsne, 'cluster').get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)
        tsne_citeseq = Tsne_dropdown_plot(df_tsne,'Citeseq',feature_name_list).get_plotly_div()
        self.add_data(tsne_citeseq=tsne_citeseq)