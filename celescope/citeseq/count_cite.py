import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import glob
import scanpy as sc
import anndata as ad

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.__init__ import HELP_DICT
from celescope.tools.matrix import CountMatrix, Features


TAG_COL = 'tag_name'


def filter_fun(df):
    """
    Filter data according to mean and standard deviation
    """
    df['filter_num'] = np.mean(df,axis=1)+np.std(df,axis=1)
    for index,_ in df.iterrows():
        df.loc[index] = df.loc[index].mask(df.loc[index]<df.loc[index]['filter_num'],0)
        #df.loc[index].apply(lambda x: 0 if x>df.loc[index]['filter_num'] else x)
    df.drop(columns='filter_num',inplace=True)
    return df



class CountMatrix2(CountMatrix):
    @classmethod
    def from_dataframe2(cls, df,features: Features, barcodes=None):
        matrix = df
        mtx = scipy.sparse.coo_matrix(matrix, matrix.shape, int)
        return cls(features,barcodes,mtx)


class Count_cite(Step):

    def __init__(self, args, display_title):
        super().__init__(args, display_title)

        self.df_read_count = pd.read_csv(args.read_count_file, sep='\t', index_col=0)

        self.match_dict = utils.parse_match_dir(args.match_dir)
        self.match_barcode = self.match_dict['match_barcode']
        self.match_matrix_dir = self.match_dict['matrix_dir']

        # input
        self.tsne_coord = glob.glob(f'{args.match_dir}/*analysis*/*tsne_coord.tsv')[0]

        # out
        self.mtx = f'{self.out_prefix}_citeseq.mtx.gz'
        self.raw_matrix_dir = f'{self.out_prefix}_raw_citeseq_matrix'
        self.filtered_matrix_dir  = f'{self.out_prefix}_filtered_citeseq_matrix'
        self.filtered_tsne_coord = f'{self.out_prefix}_filtered_tsne_coord.tsv'
        

    @utils.add_log
    def run(self):
        mapped_read = int(self.df_read_count['read_count'].sum())

        # in cell
        df_read_count_in_cell = self.df_read_count[self.df_read_count.index.isin(self.match_barcode)]
        mapped_read_in_cell = int(df_read_count_in_cell['read_count'].sum())
        self.add_metric(
            name='Mapped Reads in Cells',
            value=mapped_read_in_cell,
            total=mapped_read,
        )

        tag_names = df_read_count_in_cell[TAG_COL].unique()
        features_raw = Features(tag_names)
        
        # raw_matrix
        raw_citeseq_matrix = CountMatrix.from_dataframe(df_read_count_in_cell, features_raw, barcodes=self.match_barcode, row=TAG_COL, column='barcode', value='UMI')
        rna_matrix = CountMatrix.from_matrix_dir(matrix_dir=self.match_matrix_dir)
        raw_merged_matrix = rna_matrix.concat_by_barcodes(raw_citeseq_matrix)
        raw_merged_matrix.to_matrix_dir(self.raw_matrix_dir)

        # UMI
        df_UMI_in_cell = df_read_count_in_cell.reset_index().groupby([
            'barcode', TAG_COL]).agg({'UMI': 'count'})

        df_UMI_in_cell = df_UMI_in_cell.reset_index()
        df_UMI_in_cell = df_UMI_in_cell.pivot(
            index='barcode', columns=TAG_COL, values='UMI')
        df_cell = pd.DataFrame(index=self.match_barcode)
        df_UMI_cell = pd.merge(
            df_cell,
            df_UMI_in_cell,
            how="left",
            left_index=True,
            right_index=True)

        # fillna
        df_UMI_cell.fillna(0, inplace=True)
        df_UMI_cell = df_UMI_cell.astype(int)
        df_UMI_cell_out = df_UMI_cell.T
        df_UMI_cell_out.to_csv(self.mtx, sep='\t', compression='gzip')


        # filter_matrix
        df_filtered = filter_fun(df_UMI_cell_out)
        features_filtered = Features(df_filtered.index.to_list())
        filtered_citeseq_matrix = CountMatrix2.from_dataframe2(df_filtered, features_filtered, barcodes=self.match_barcode)
        filtered_merged_matrix = rna_matrix.concat_by_barcodes(filtered_citeseq_matrix)
        filtered_merged_matrix.to_matrix_dir(self.filtered_matrix_dir)


        # normalize citeseq date 
        obs = pd.DataFrame(index=df_filtered.columns)
        var = pd.DataFrame(df_filtered.index,index=df_filtered.index,columns=['ADT'])
        mdata_citeseq = ad.AnnData(np.array(df_filtered.T),obs=obs,var=var)
        sc.pp.normalize_total(
                mdata_citeseq,
                target_sum=1e4,
                inplace=True,
            )
        sc.pp.log1p(mdata_citeseq)

        # filtered_tsne.csv
        df_tsne = pd.read_csv(self.tsne_coord,sep="\t")
        if 'Unnamed: 0' in df_tsne.columns:
            df_tsne.rename(columns={'Unnamed: 0': 'barcode'}, inplace=True)
            df_tsne = df_tsne.set_index('barcode')

        df_citeseq = mdata_citeseq.to_df()
        df_tsne = pd.concat([df_tsne,df_citeseq],axis=1)
        df_tsne.fillna(0,inplace=True)
        df_tsne.to_csv(self.filtered_tsne_coord,sep='\t')

        

        # UMI
        UMIs = df_UMI_cell.apply(sum, axis=1)
        median_umi = round(np.median(UMIs), 2)
        mean_umi = round(np.mean(UMIs), 2)
        self.add_metric(
            name='Median UMI per Cell',
            value=float(median_umi),
        )
        self.add_metric(
            name='Mean UMI per Cell',
            value=float(mean_umi),
        )


def get_opts_count_cite(parser, sub_program):
    if sub_program:
        parser.add_argument("--match_dir", help=HELP_DICT['match_dir'], required=True)
        parser.add_argument("--read_count_file", help="tag read count file")
        s_common(parser)


def count_cite(args):

    with Count_cite(args, display_title='Count') as runner:
        runner.run()
