import unittest
import os
import pandas as pd
from celescope.snp.snpCalling import convert
from celescope.tools.utils import read_barcode_file


class test_snp(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/snp')
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92'
        self.gene_list_file = './gene_list.tsv'

    def test_convert(self):
        gene_list = convert(self.gene_list_file, self.genomeDir)
        print(gene_list)

if __name__ == '__main__':
    unittest.main()
