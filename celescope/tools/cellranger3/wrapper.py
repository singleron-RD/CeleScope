"""
cellranger3 cell-calling on raw matrix
"""

import os
import sys

import numpy as np
import pandas as pd
import scipy.io

from celescope.tools.__init__ import (BARCODE_FILE_NAME, FEATURE_FILE_NAME,
                                      MATRIX_FILE_NAME)
from celescope.tools.cellranger3.cell_calling_3 import find_nonambient_barcodes
import celescope.tools.utils as utils


EXPECTED_CELL_NUM = 3000


@utils.add_log
def read_raw_matrix(all_matrix_10X_dir):

    raw_mat_path = os.path.join(all_matrix_10X_dir, MATRIX_FILE_NAME)
    raw_mat = scipy.io.mmread(raw_mat_path)  # scipy.sparse.coo.coo_matrix

    raw_features_path = os.path.join(all_matrix_10X_dir, FEATURE_FILE_NAME)

    raw_barcodes_path = os.path.join(all_matrix_10X_dir, BARCODE_FILE_NAME)
    raw_barcodes_df = pd.read_csv(raw_barcodes_path, sep='\t', error_bad_lines=False, names=['barcode'])
    raw_barcodes = np.array(raw_barcodes_df['barcode'].tolist())

    return raw_mat, raw_features_path, raw_barcodes


@utils.add_log
def cell_calling(raw_mat, expected_cell_num):

    filtered_bc_indices, _round_1_filtered_metrics, _non_ambient_barcode_result = find_nonambient_barcodes(
        raw_mat=raw_mat, recovered_cells=expected_cell_num)

    return filtered_bc_indices


class Cell_calling():
    def __init__(self, outdir, raw_mat, raw_features_path, raw_barcodes):
        self.raw_mat = raw_mat
        self.raw_features_path = raw_features_path
        self.raw_barcodes = raw_barcodes
        self.expected_cell_num = EXPECTED_CELL_NUM
        self.out_mat = f'{outdir}/{MATRIX_FILE_NAME}'
        self.out_feature = f'{outdir}/{FEATURE_FILE_NAME}'
        self.out_barcode = f'{outdir}/{BARCODE_FILE_NAME}'
        utils.check_mkdir(outdir)

    def write_slice_matrix(self, filtered_bc_indices):
        mtx_csc = self.raw_mat.tocsc()
        sliced_mtx = mtx_csc[:, filtered_bc_indices]
        scipy.io.mmwrite(target=self.out_mat, a=sliced_mtx)

        barcodes = pd.Series(self.raw_barcodes[filtered_bc_indices])
        barcodes.to_csv(self.out_barcode, index=False, sep='\t', header=False)

        os.system(f"cp {self.raw_features_path} {self.out_feature}")

    @utils.add_log
    def run(self):
        filtered_bc_indices = cell_calling(self.raw_mat, self.expected_cell_num)
        self.write_slice_matrix(filtered_bc_indices)


def main():
    outdir = sys.argv[1]
    all_matrix_10X_dir = sys.argv[2]
    raw_mat, raw_features_path, raw_barcodes = read_raw_matrix(all_matrix_10X_dir)
    runner = Cell_calling(outdir, raw_mat, raw_features_path, raw_barcodes)
    runner.run()


if __name__ == '__main__':
    main()
