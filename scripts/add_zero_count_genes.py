"""
Add zeor count genes to matrix
"""
import argparse
import os

import scipy

from celescope.tools.matrix import CountMatrix, Features
from celescope.tools.reference import GtfParser
from celescope.tools import utils


def add_zero_count(cm: CountMatrix, all_features: Features):
    """
    Returns:
        CountMatrix object
    """
    m = cm.get_matrix()

    all_gene_id = all_features.gene_id
    feature_index_dict = {} # {'ENSG00000223972': 0, 'ENSG00000227232': 1,}
    for index, gene_id in enumerate(all_gene_id):
        feature_index_dict[gene_id] = index

    row_index_dict = {} # {0: 25969, 1: 54488}
    for sub_index, gene_id in enumerate(cm.get_features().gene_id):
        row_index_dict[sub_index] = feature_index_dict[gene_id]

    new_row = [row_index_dict[sub_index] for sub_index in m.row]
    new_matrix = scipy.sparse.coo_matrix((m.data, (new_row, m.col)), 
            shape=(len(feature_index_dict), len(cm.get_barcodes())))

    return CountMatrix(all_features, cm.get_barcodes(), new_matrix)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add zeor count genes to matrix')
    parser.add_argument('--matrix_list_file', help='Plain text file with one matrix dir per line.', required=True)
    parser.add_argument('--gtf', required=True, help='GTF file.')
    parser.add_argument('--outdir', help='Output directory.', default='./out')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    gp = GtfParser(args.gtf)
    gp.get_id_name()
    all_features = gp.get_features()

    matrix_list, _ = utils.read_one_col(args.matrix_list_file)
    for matrix_dir in matrix_list:
        matrix_name = matrix_dir.split('/')[-1]
        out_matrix_dir = f'{args.outdir}/{matrix_name}'
        cm = CountMatrix.from_matrix_dir(matrix_dir)
        new_cm = add_zero_count(cm, all_features)
        new_cm.to_matrix_dir(out_matrix_dir)
