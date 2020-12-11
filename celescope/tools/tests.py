import unittest
import os
import pandas as pd
from celescope.tools.STAR import Step_mapping
from celescope.tools.utils import read_barcode_file, gene_convert, get_fq
from .Chemistry import Chemistry
from celescope.tools.count import *

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
            results.append(ch.check_chemistry())
        print(results)
        assert results == ['scopeV2.1.1', 'scopeV2.1.1', 'scopeV2.0.1', 'scopeV2.2.1']
        fq = '/SGRNJ03/DATA03/2004/20201207_10/R2011312-PCRD-201203-1_combined_R1.fastq.gz'
        ch = Chemistry(fq)
        print(ch.check_chemistry())

    @unittest.skip('pass')
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

    @unittest.skip('pass')
    def test_matrix_10X(self):
        match_dir = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/test2/'
        validated_barcodes, _ncell = read_barcode_file(match_dir)    
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/')
        df = pd.read_csv('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/test2/05.count/test2_count_detail.txt.gz',sep='\t')
        outdir = 'test2/05.count'
        sample = 'test2'
        gtf_file = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf'
        matrix_10X(df, outdir, sample, gtf_file, dir_name='matrix_10X_new', validated_barcodes=validated_barcodes)

    @unittest.skip('pass')
    def test_call_cells(self):
        df = pd.read_csv('/SGRNJ02/RandD4/RD2019016/20201209/J-Demo_Y1/05.count/J-Demo_Y1_count_detail.txt.gz',sep='\t')
        os.chdir = '/SGRNJ02/RandD4/RD2019016/20201209/'
        sample = 'J-Demo_Y1'
        outdir = f'{sample}/05.count'
        cells = 'auto'
        pdf =  outdir + '/barcode_filter_magnitude.pdf'
        marked_counts_file = outdir + '/' + sample + '_counts.txt'
        (validated_barcodes, threshold, cell_num, CB_describe) = call_cells(
            df, cells, pdf, marked_counts_file)
        print(threshold)

    @unittest.skip('pass')
    def test_rescue(self):
        dir_name = 'all_matrix'
        os.chdir('/SGRNJ02/RandD4/RD2019016/20201209/')
        sample = 'J-Demo_Y1'
        outdir = f'{sample}/05.count'
        threshold = 279
        matrix_dir = f"{outdir}/{sample}_{dir_name}/"
        threshold = rescue_cells(outdir, sample, matrix_dir, threshold)
        print(threshold)
        df = pd.read_csv('/SGRNJ02/RandD4/RD2019016/20201209/J-Demo_Y1/05.count/J-Demo_Y1_count_detail.txt.gz',sep='\t')
        df_sum = df.groupby(['Barcode']).agg({'UMI':'count'})
        validated_barcodes = get_validated_barcodes(df_sum, threshold, col='UMI')
        print(len(validated_barcodes))

    def test_count_pipe(self):
        count_detail_file = '/SGRNJ02/RandD4/RD2019016/20201209/J-Demo_Y1/05.count/J-Demo_Y1_count_detail.txt.gz'
        df = pd.read_csv(count_detail_file, sep='\t')
        dir_name = 'all_matrix'
        os.chdir('/SGRNJ02/RandD4/RD2019016/20201209/')
        sample = 'J-Demo_Y1'
        outdir = f'{sample}/05.count'
        rescue = True
        gtf_file = '/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.92.chr.gtf'
        cells = 'auto'
        assay = 'rna'

        # call cells
        pdf = outdir + '/barcode_filter_magnitude.pdf'
        df_sum, threshold = call_cells(df, cells, pdf)

        # rescue low UMI cells
        if rescue:
            matrix_dir = f"{outdir}/{sample}_{dir_name}/"
            threshold = rescue_cells(outdir, sample, matrix_dir, threshold)
        
        # get cell stats
        marked_counts_file = outdir + '/' + sample + '_counts.txt'
        validated_barcodes, CB_describe = get_cell_stats(df_sum, threshold, marked_counts_file)

        # export cell matrix
        matrix_10X(df, outdir, sample, gtf_file, dir_name='matrix_10X', validated_barcodes=validated_barcodes)
        (CB_total_Genes, CB_reads_count,
            reads_mapped_to_transcriptome) = expression_matrix(
                df, validated_barcodes, outdir, sample, gtf_file)

        # downsampling
        validated_barcodes = set(validated_barcodes)
        downsample_file = outdir + '/' + sample + '_downsample.txt'
        Saturation = downsample(
            count_detail_file,
            validated_barcodes,
            downsample_file)

        # summary
        stat_file = outdir + '/stat.txt'
        get_summary(df, sample, Saturation, CB_describe, CB_total_Genes,
                    CB_reads_count, reads_mapped_to_transcriptome, stat_file,
                    outdir + '/../')

        report_prepare(marked_counts_file, downsample_file, outdir + '/..')

        t = reporter(assay=assay,
                    name='count', sample=sample,
                    stat_file=outdir + '/stat.txt',
                    outdir=outdir + '/..')
        t.get_report()

if __name__ == '__main__':
    unittest.main()