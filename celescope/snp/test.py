import unittest
import os
import pandas as pd
import glob
from celescope.snp.snpCalling import *
from celescope.tools.utils import read_barcode_file, glob_genomeDir, parse_vcf
from .analysis_snp import analysis_variant


class test_snp(unittest.TestCase):
    '''
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/snp')
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92'
        self.gene_list_file = './gene_list.tsv'
        self.index_file = './S20070818_TS/05.snpCalling/S20070818_TS_cell_index.tsv'
        self.thread = 20
        self.outdir = './S20070818_TS/05.snpCalling/'
        self.analysis_outdir = './S20070818_TS/06.analysis_snp'
        _refFlat, self.gtf, self.fasta = glob_genomeDir(self.genomeDir, fa=True)
        self.match_dir = '/SGRNJ02/RandD4/RD20051303_Panel/20200717/S20070818_ZL/'
        self.sample = 'S20070818_TS'
        self.count_file = './S20070818_TS/05.snpCalling/S20070818_TS_count.tsv'
        self.vcf_file = './S20070818_TS/05.snpCalling/S20070818_TS_anno.vcf'
        self.assay = 'snp'
        self.annovar_config = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/soft/annovar/annovar.config'
        self.step_analysis_variant = analysis_variant(
            self.analysis_outdir,
            self.sample,
            self.match_dir,
            self.vcf_file,
            self.index_file,
            self.assay,
            self.annovar_config,
        )
    '''

    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/snp')
        sample = '20201224_consplit'
        self.sample = sample
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92'
        self.gene_list_file = './gene_list.tsv'
        self.CID_file = f'./{sample}/05.snpCalling/{sample}_CID.tsv'
        self.thread = 20
        self.analysis_outdir = f'./{sample}/06.analysis_snp'
        _refFlat, self.gtf, self.fasta = glob_genomeDir(self.genomeDir, fa=True)
        self.match_dir = '/SGRNJ02/RandD4/RD20051303_Panel/20200717/S20070818_ZL/'
        self.vcf_file = f'./{sample}/05.snpCalling/{sample}_merged.vcf'
        self.assay = 'snp'
        self.annovar_config = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/soft/annovar/annovar.config'
        self.step_analysis_variant = analysis_variant(
            self.analysis_outdir,
            self.sample,
            self.match_dir,
            self.vcf_file,
            self.CID_file,
            self.assay,
            self.annovar_config,
        )


    @unittest.skip('pass')
    def test_annovar(self):
        self.step_analysis_variant.annovar()
    
    @unittest.skip('pass')
    def test_parse_annovar(self):
        from celescope.tools.utils import parse_annovar
        self.annovar_file = glob.glob(f'{self.analysis_outdir}/{self.sample}*_multianno.txt')[0]
        parse_annovar(self.annovar_file)

    def test_parse_vcf(self):
        self.step_analysis_variant.get_df_vcf()
        print(self.step_analysis_variant.df_vcf)
        

if __name__ == '__main__':
    unittest.main()
