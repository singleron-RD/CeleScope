import scipy.io
import scipy.sparse
import pandas as pd

from celescope.tools.__init__ import (BARCODE_FILE_NAME, FEATURE_FILE_NAME, MATRIX_FILE_NAME)
from celescope.tools import utils

class Features:
    def __init__(self, gene_id: list, gene_name=None, type=None):
        """
        Args:
            gene_id: list of gene id
            gene_name: list of gene name
            type: ype of features, e.g. [gene, protein]
        """
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.type = type

    @classmethod
    def from_tsv(cls, tsv_file):
        df = pd.read_csv(tsv_file, sep='\t', on_bad_lines='skip', names=['gene_id', 'gene_name', 'type'])
        gene_id = df['gene_id'].tolist()
        gene_name = df['gene_name'].tolist()
        type = df['type'].tolist()
        return cls(gene_id, gene_name, type)

    def to_tsv(self, tsv_file):
        df = pd.DataFrame({'gene_id': self.gene_id, 'gene_name': self.gene_name, 'type': self.type})
        df.to_csv(tsv_file, sep='\t', index=False, header=False)


class CountMatrix:
    def __init__(self, features: Features, barcodes: list, matrix):
        """
        Args:
            features: Features object
            barcodes: list of barcodes
            matrix: scipy.sparse.coo_matrix
        """
        self.__features = features
        self.__barcodes = barcodes
        self.__matrix = matrix
        self.shape = matrix.shape

    @classmethod
    @utils.add_log
    def from_matrix_dir(cls, matrix_dir):
        features_tsv = f'{matrix_dir}/{FEATURE_FILE_NAME}'
        features = Features.from_tsv(tsv_file=features_tsv)
        barcode_file = f'{matrix_dir}/{BARCODE_FILE_NAME}'
        barcodes, _ = utils.read_one_col(barcode_file)
        matrix_path = f'{matrix_dir}/{MATRIX_FILE_NAME}'
        matrix= scipy.io.mmread(matrix_path)

        return cls(features, barcodes, matrix)

    def to_matrix_dir(self, matrix_dir):
        utils.check_mkdir(dir_name=matrix_dir)
        features = self.__features.to_tsv(f'{matrix_dir}/{FEATURE_FILE_NAME}')
        barcodes = pd.Series(self.__barcodes).to_csv(f'{matrix_dir}/{BARCODE_FILE_NAME}', index=False, sep='\t', header=False)
        matrix_path = f'{matrix_dir}/{MATRIX_FILE_NAME}'
        scipy.io.mmwrite(matrix_path, self.__matrix)

    @classmethod
    def from_dataframe(cls, df, row='geneID', column='Barcode', value="UMI", gtf_dict=None, type=None):
        """
        Args:
            df: dataframe with columns: [row, column, value]. Will be grouped by row and column and count lines of value.
            value: value name in df, UMI
            gtf_dict: {gene_id: gene_name}
            type: type of features, e.g. [gene, protein]
        """
        df = df.groupby([row, column]).agg({value: 'count'})
        mtx = scipy.sparse.coo_matrix((df[value], (df.index.codes[0], df.index.codes[1])))
        gene_id = df.index.levels[0].tolist()
        # add gene symbol
        gene_name = None
        if gtf_dict:
            gene_name = [gtf_dict[x] for x in gene_id]

        features = Features(gene_id, gene_name, type)
        barcodes = df.index.levels[1].tolist()
        return cls(features, barcodes, mtx)

    def __str__(self):
        n_row, n_col = self.shape[0], self.shape[1]
        return f"CountMatrix object\n {n_row} x {n_col} coo_matrix"


    def __repr__(self):
        return self.__str__()

    def concat_by_barcodes(self, other):
        if other.get_barcodes() != self.get_barcodes():
            raise ValueError('barcodes are not the same')

        if inter := set(self.get_features().gene_id).intersection(set(other.get_features().gene_id)):
            raise ValueError(f'deuplicated gene_id: {inter}')

        gene_id = self.get_features().gene_id + other.get_features().gene_id
        gene_name = self.get_features().gene_name + other.get_features().gene_name
        type = None
        if self.get_features().type and other.get_features().type:
            type = self.get_features().type + other.get_features().type
        features = Features(gene_id, gene_name, type)

        matrix = scipy.sparse.vstack([self.get_matrix(), other.get_matrix()])

        return CountMatrix(features, self.__barcodes, matrix)

    def get_barcodes(self):
        return self.__barcodes

    def get_features(self):
        return self.__features
    
    def get_matrix(self):
        return self.__matrix

