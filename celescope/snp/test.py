import unittest
import os
import pandas as pd
from celescope.snp.snpCalling import convert, call_all_snp
from celescope.tools.utils import read_barcode_file, glob_genomeDir


class test_snp(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/snp')
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92'
        self.gene_list_file = './gene_list.tsv'
        self.index_file = './S20070818_TS/05.snpCalling/S20070818_TS_cell_index.tsv'
        self.thread = 4
        self.outdir = './S20070818_TS/05.snpCalling/'

    @unittest.skip('pass')
    def test_convert(self):
        _refFlat, gtf = glob_genomeDir(self.genomeDir)
        gene_list = convert(self.gene_list_file, gtf)
        print(gene_list)

    def test_call_all_snp(self):
        _refFlat, gtf, fasta = glob_genomeDir(self.genomeDir, fa=True)
        call_all_snp(self.index_file, self.outdir, self.thread, fasta)


if __name__ == '__main__':
    unittest.main()
