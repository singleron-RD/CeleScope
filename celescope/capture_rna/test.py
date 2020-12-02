import unittest
import os
import pandas as pd
from celescope.capture_rna.count_capture_rna import expression_matrix
from celescope.tools.utils import read_barcode_file, read_one_col, glob_genomeDir
from celescope.tools.report import reporter


class testHLA(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/0910_panel/')
        self.sample = 'S20071508_D_TS'
        count_detail_file = './/S20071508_D_TS/05.count_capture_rna/S20071508_D_TS_count_detail.txt'
        self.df = pd.read_table(count_detail_file, header=0)
        self.match_dir = '/SGRNJ02/RandD4/RD20051303_Panel/20200729/S20071508_D_ZL'
        self.sc_cell_barcodes, self.sc_cell_number = read_barcode_file(self.match_dir)
        self.outdir = f'{self.sample}/05.count_capture_rna/'
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92'
        self.validated_barcodes, _ = read_one_col(f'{self.sample}/05.count_capture_rna/{self.sample}_matrix_10X/barcodes.tsv') 
        _refFlat, self.gtf = glob_genomeDir(self.genomeDir)
        self.assay = 'capture_rna'

    #@unittest.skip('pass')
    def test_expression_matrix(self):
        (
            CB_total_Genes,
            CB_reads_count,
            reads_mapped_to_transcriptome,
            match_cell_str,
            match_UMI_median
        ) = expression_matrix(
            self.df,
            self.validated_barcodes,
            self.outdir,
            self.sample,
            self.gtf,
            self.sc_cell_barcodes,
            self.sc_cell_number
            )
        print(match_cell_str)
        print(match_UMI_median)

    def test_report(self):
        t = reporter(assay=self.assay,
            name='count_capture_rna', sample=self.sample,
            stat_file=self.outdir + '/stat.txt',
            outdir=self.outdir + '/..')
        t.get_report()




if __name__ == '__main__':
    unittest.main()