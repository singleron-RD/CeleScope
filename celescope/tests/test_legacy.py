import unittest
import os
import pandas as pd
from collections import namedtuple
from celescope.tools.STAR import Step_mapping
from celescope.tools.utils import *
from .Chemistry import Chemistry
from celescope.tools.count import *
from celescope.tools.count import report_prepare as report_p
from celescope.tools.cutadapt import *
from celescope.tools.Multi import *

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
        os.chdir('/SGRNJ03/PROJ03/PROJ_20.SC/PN20122407_SCOPEv2/temp/CK2101_origin/')
        sample = 'CK2101'
        outdir = f'{sample}/05.count'
        df = pd.read_csv(f'{outdir}/{sample}_count_detail.txt.gz', sep='\t')
        outdir = f'{sample}/05.count'
        cell_calling_method = 'cellranger3'
        force_cell_num = None
        expected_cell_num = 3000
        all_matrix_10X_dir = f'{outdir}/{sample}_all_matrix'

        # df_sum
        df_sum = get_df_sum(df)

        # call cells
        cell_bc, threshold = cell_calling(
            cell_calling_method, force_cell_num, expected_cell_num, all_matrix_10X_dir, df_sum, outdir, sample)

        # plot
        cell_num = len(cell_bc)
        print(f'cell_num:{cell_num}')
        #plot_barcode_UMI(df_sum, threshold, expected_cell_num, cell_num, outdir, sample, cell_calling_method,col='UMI')


    @unittest.skip('pass')
    def test_read_adapter_fasta(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna')
        adapter_fasta = './adapter.fasta'
        adapter_args = read_adapter_fasta(adapter_fasta)
        assert adapter_args == ['a1=ATCG','a2=TGCAA']

    @unittest.skip('pass')
    def test_chemistry_jiace(self):
        fq = '/SGRNJ03/DATA03/2004/20201224_18/R20067441-PCRD-201207-1-R2012093_combined_R1.fastq.gz\
,/SGRNJ03/DATA03/2004/20201228_6/R2012093-PCRD-201207-1_combined_R1.fastq.gz'
        ch = Chemistry(fq)
        print(ch.check_chemistry())

    '''
    @unittest.skip('pass')
    def test_report_prepare(self):
        outdir = "/SGRNJ03/randd/P19112803_SCOPEv1/test1/NJXK01_1/05.count"
        count_file = f"{outdir}/NJXK01_1_counts.txt"
        downsample_file = f"{outdir}/NJXK01_1_downsample.txt"
        report_prepare(count_file, downsample_file, outdir)
    '''

    @unittest.skip('pass')
    def test_Multi(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/rebuild')
        multi = Multi('rna')

        sys.argv = ['celescope', '--mapfile', 'test.mapfile', '--not_gzip']
        '''
        multi.parse_args()
        print(multi.parse_step_args('sample'))
        multi.prepare()
        multi.sample('test1')
        multi.barcode('test1')
        print(multi.sjm_cmd)
        '''

        multi.run()
        print(multi.sjm_cmd)

    @unittest.skip('pass')
    @add_mem
    def test_downsample(self):
        assay = 'rna'
        step = 'count'
        outdir = '/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/rebuild/test1/05.count/'
        sample = 'test1'
        count_detail_file = f'{outdir}/{sample}_count_detail.txt.gz'
        cell_bc, _ = read_barcode_file(f"{outdir}/../")
        df = pd.read_table(count_detail_file, header=0)
        df_cell = df.loc[df["Barcode"].isin(cell_bc), :]

        json_file = f'{outdir}/data.json'
        reporter = Reporter(assay, 'count', sample, outdir, json_file)
        downsample_file = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/rebuild/test1/05.count/new.downsample"
        Saturation, res_dict = downsample(df, cell_bc, downsample_file)
        print(res_dict)
        reporter.add_data_item(downsample_summary=res_dict)

        reporter.dump_json()
    
    @unittest.skip('pass')
    @add_mem
    def test_downsample_large(self):
        count_detail_file = '/SGRNJ03/randd/P19112803_SCOPEv1/test1/NJXK01_1/old.count/NJXK01_1_count_detail.txt.gz'
        cell_bc, _ = read_barcode_file("/SGRNJ03/randd/P19112803_SCOPEv1/test1/NJXK01_1/")
        df = pd.read_table(count_detail_file, header=0)
        #df_cell = df.loc[df["Barcode"].isin(cell_bc), :]

        downsample_file = "/SGRNJ03/randd/P19112803_SCOPEv1/test1/NJXK01_1/05.count/test.downsample"
        downsample(df, cell_bc, downsample_file)
    
    @unittest.skip('pass')
    @add_mem
    def test_downsample_real1(self):
        count_detail_file = '/SGRNJ03/PROJ03/PROJ_20.SC/PN20122407_SCOPEv2/temp/CK2101_origin/CK2101/05.count/CK2101_count_detail.txt.gz'
        cell_bc, _ = read_barcode_file("/SGRNJ03/PROJ03/PROJ_20.SC/PN20122407_SCOPEv2/temp/CK2101_origin/CK2101/")
        df = pd.read_table(count_detail_file, header=0)
        #df_cell = df.loc[df["Barcode"].isin(cell_bc), :]

        downsample_file = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/downsample/CK2101_origin.downsample.txt"
        downsample(df, cell_bc, downsample_file)

    @unittest.skip('pass')
    @add_mem
    def test_downsample_large1(self):
        assay = 'rna'
        step = 'count'
        outdir = '/SGRNJ03/PROJ03/PROJ_20.SC/PN20122407_SCOPEv2/temp/CK2101_origin/CK2101/05.count/'
        sample = 'CK2101'
        count_detail_file = f'{outdir}/{sample}_count_detail.txt.gz'
        cell_bc, _ = read_barcode_file(f"{outdir}/../")
        df = pd.read_table(count_detail_file, header=0)
        df_cell = df.loc[df["Barcode"].isin(cell_bc), :]

        json_file = f'/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/downsample/CK2101_origin.json'
        reporter = Reporter(assay, 'count', sample, outdir, json_file)
        downsample_file = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/downsample/CK2101_origin.downsample.txt"
        Saturation, res_dict = downsample(df, cell_bc, downsample_file)
        print(res_dict)
        reporter.add_data_item(downsample_summary=res_dict)

        reporter.dump_json()

    @unittest.skip('pass')
    @add_mem
    def test_downsample_large2(self):
        assay = 'rna'
        step = 'count'
        outdir = '/SGRNJ03/PROJ03/PROJ_20.SC/PN20122407_SCOPEv2/temp/CK816_1-1/CK816_1-1/05.count/'
        sample = 'CK816_1-1'
        count_detail_file = f'{outdir}/{sample}_count_detail.txt.gz'
        cell_bc, _ = read_barcode_file(f"{outdir}/../")
        df = pd.read_table(count_detail_file, header=0)
        df_cell = df.loc[df["Barcode"].isin(cell_bc), :]

        json_file = f'/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/downsample/CK816_1-1.json'
        reporter = Reporter(assay, 'count', sample, outdir, json_file)
        downsample_file = "/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/downsample/CK816_1-1.downsample.txt"
        Saturation, res_dict = downsample(df, cell_bc, downsample_file)
        print(res_dict)
        reporter.add_data_item(downsample_summary=res_dict)

        reporter.dump_json()

    @unittest.skip('pass')
    @add_mem
    def test_report_prepare(self):
        os.chdir("/SGRNJ03/randd/P19112803_SCOPEv1/test1/")
        marked_counts_file = 'NJXK01_1/05.count/NJXK01_1_counts.txt'
        downsample_file = 'NJXK01_1/05.count/NJXK01_1_downsample.txt'
        outdir = 'NJXK01_1/05.count/'
        report_p(marked_counts_file, downsample_file, outdir)

    @unittest.skip('pass')
    def test_report(self):
        os.chdir("/SGRNJ03/randd/P19112803_SCOPEv1/test1/")
        assay = 'rna'
        sample = 'NJXK01_1'
        outdir = 'NJXK01_1/05.count/'
        t = reporter(assay=assay,
                 name='count', sample=sample,
                 stat_file=outdir + '/stat.txt',
                 outdir=outdir + '/..')
        t.get_report()

    @unittest.skip('pass')
    def test_metrics(self):
        assay = 'rna'
        sample = 'test1'
        __STEPS__ = [
            'sample',
            'barcode',
            'cutadapt',
            'STAR',
            "featureCounts",
            "count",
            'analysis',
        ]
        for index, step in enumerate(__STEPS__):
            print(index, step)
            outdir = f'/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/rna/rebuild/test1/0{index}.{step}'
            report = Reporter(
                assay,
                step,
                sample,
                outdir,
            )
            report.stat_to_json()
            report.dump_json()

    def test_sorted_dumb_consensus(self):
        os.chdir('/SGRNJ01/RD_dir/pipeline_test/zhouyiqi/unittest/vdj/rebuild/TCR_20201221_consensus/03.consensus')
        fastq = 'TCR_20201221_consensus_sorted.fq.tmp'
        sorted_dumb_consensus(fastq, outfile = 'update_consensus.fq.gz', threshold=0.5)