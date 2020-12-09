import unittest
import os
import pandas as pd
from celescope.tools.STAR import Step_mapping
from celescope.tools.utils import read_barcode_file, gene_convert, get_fq
from .Chemistry import Chemistry
import subprocess

class Tests(unittest.TestCase):
    def setUp(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna')
        self.sample = 'test1'
        self.fq = './test1/02.cutadapt/test1_clean_2.fq.gz'
        self.genomeDir = '/SGRNJ/Public/Database/genome/homo_mus'
        self.thread = 4
        self.assay = 'rna'
        self.out_unmapped = False
        self.debug = True

    @unittest.skip('pass')
    def test_ribo(self):
        self.outdir = f'{self.sample}/03.STAR'
        mapping = Step_mapping(
            self.sample, 
            self.outdir, 
            self.assay, 
            self.thread,
            self.fq, 
            self.genomeDir, 
            self.out_unmapped, 
            self.debug)
        mapping.STAR_map_log = f'{self.outdir}/{self.sample}_Log.final.out'  
        mapping.STAR_map_log = f'{self.outdir}/{self.sample}_Log.final.out'
        mapping.picard_region_log = f'{self.outdir}/{self.sample}_region.log'
        mapping.ribo()
        mapping.format_stat()
        mapping.report()

    @unittest.skip('pass')
    def test_chemistry(self):
        fq1s = [
            '/SGRNJ/DATA_PROJ/2004/20201102_9/MG_201026_1-R2010283_L3_1.fq.gz',
            '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/test_data/smk/R2009217_L4_1.fq.gz',
            '/SGRNJ/DATA_PROJ/2003/20200727_6/R2007185_L3_1.fq.gz',
            '/SGRNJ/DATA_PROJ/2004/20201029/MG_201023_1-R2010226_L3_1.fq.gz',
        ]
        results = []
        for fq1 in fq1s:
            ch = Chemistry(fq1)
            results.append(ch.get_chemistry())
        print(results)
        assert results == ['scopeV2.1.1', 'scopeV2.1.1', 'scopeV2.0.1', 'scopeV2.2.1']
        fq = '/SGRNJ03/DATA03/2004/20201122_14/R2011312_R1.fastq.gz'
        ch = Chemistry(fq)
        print(ch.get_chemistry())

    #@unittest.skip('pass')
    def test_gtf(self):
        '''
        gtf_file = '/SGRNJ/Database/script/genome/hs/gtf/Homo_sapiens.GRCh38.99.gtf'
        id_name = gene_convert(gtf_file)
        print(f"ENSG00000001629: {id_name['ENSG00000001629']}")

        gtf_file = '/SGRNJ01/RD_dir/pipeline_test/litao/genomes/Cricetulus_griseus/Cricetulus_griseus_crigri.CriGri_1.0.101.gtf'
        id_name = gene_convert(gtf_file)
        print(id_name)
        '''

        gtf_file = '/Public/Database/genome/Sus_scrofa/ncbi/GCF_000003025.6_Sscrofa11.1_genomic_new.gtf'
        id_name = gene_convert(gtf_file)
        print(id_name['tRNA-Asp'])

    @unittest.skip('pass')
    def test_rescue(self):
        tools_dir = os.path.dirname(__file__)
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna')
        outdir = 'test2/05.count'
        sample = 'test2'
        matrix_dir = f'{outdir}/{sample}_all_matrix_10X'
        app = f'{tools_dir}/rescue.R'
        cmd = (
            f'Rscript {app} '
            f'--outdir {outdir} '
            f'--sample {sample} '
            f'--matrix_dir {matrix_dir} '
            f'--threshold 14 '
        )
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    
    @unittest.skip('pass')
    def test_get_fq(self):
        library_id = 'R2011312'
        library_path = '/SGRNJ03/DATA03/2004/20201122_14/'
        fq1, fq2 = get_fq(library_id, library_path)
        print(fq1, fq2)


if __name__ == '__main__':
    unittest.main()