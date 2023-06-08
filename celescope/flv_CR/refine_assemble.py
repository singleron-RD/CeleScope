from xopen import xopen
from celescope.tools import utils

import os
import subprocess
import pysam
import pandas as pd


class Filter_noise:
    """ Filter noise barcode.
    For is_cell=False barcodes, keep the highest UMI chain when its UMI >= coeff * second largest UMI. Then only paired-chain barcodes are keeped.

    First sort in descending order according to barcode, chain, umis, so that the next step does not need to be sorted Perform groupby by barcode and chain.
    If each groupby object has only one data, keep it; if it is more than 1, check whether the highest umi is greater than the second highest umi*1.5, if it is satisfied, keep the highest one, if not, then all not reserved Perform groupby by barcode. 
    If the chain number of each barcode is not equal to 2, it will not be retained.
    When retaining, use a set to record the contig_id, and after all data checks are completed, determine which rows need to be retained according to the contig_id.
    """
    
    def __init__(self, df, coeff, seqtype):

        self.df = df.sort_values("umis", ascending=False)
        self.coeff = coeff
        self.seqtype = seqtype
    
    @utils.add_log
    def __call__(self):
        
        df_sorted = self.df.sort_values(['barcode', 'chain', 'umis'], ascending=[False, False, False])

        # 按barcode和chain进行groupby，保留满足条件的数据
        grouped = df_sorted.groupby(['barcode', 'chain'])
        filtered_data = pd.DataFrame(columns=df_sorted.columns)  # 创建一个空的DataFrame来存储过滤后的数据
        contig_id_set = set()  # 用于记录需要保留的contig_id

        for _group_name, group_data in grouped:
            if len(group_data) == 1:
                contig_id_set.add(group_data['contig_id'].values[0])  
            else:
                sorted_umis = group_data['umis'].values
                highest_umi = sorted_umis[0]
                second_highest_umi = sorted_umis[1]
                if highest_umi >= second_highest_umi * self.coeff:
                    contig_id_set.add(group_data.iloc[0]['contig_id'])  # 记录保留的contig_id

        filtered_data = self.df[self.df['contig_id'].isin(contig_id_set)]

        grouped_barcode = filtered_data.groupby('barcode')
        final_filtered_data = pd.DataFrame(columns=filtered_data.columns)
        contig_id_set = set()

        for _group_name, group_data in grouped_barcode:
            if len(group_data['chain'].unique()) == 2:
                if self.seqtype == 'BCR' and 'IGH' not in group_data['chain'].unique():
                    continue
                contig_id_set.update(set(group_data['contig_id']))  # 记录保留的contig_id

        # 根据contig_id确定需要保留的行
        final_filtered_data = filtered_data[filtered_data['contig_id'].isin(contig_id_set)]
        final_filtered_data = final_filtered_data.sort_values(['barcode'])

        return final_filtered_data


class Refine_assemble:
    """
    ## Features

    - TCR/BCR Assemble by Cellranger.

    - Generate Mapping, Cells, V(D)J annotations metrics in html.

    ## Output
    
    - `03.assemble/{sample}_refine` Cellranger vdj refined results.

    """
    def __init__(self, args):

        # common
        self.sample = args.sample
        self.outdir = args.outdir
        self.coeff = args.coeff
        self.seqtype = args.seqtype
        
        # in
        self.raw_all_contig_anno = pd.read_csv(f"{self.outdir}/{self.sample}/outs/all_contig_annotations.csv")
        self.raw_filter_contig_anno = pd.read_csv(f"{self.outdir}/{self.sample}/outs/filtered_contig_annotations.csv")
        self.all_bam = f"{self.outdir}/{self.sample}/outs/all_contig.bam"
        self.fqs_dir = os.path.abspath(args.fqs_dir)
        self.raw_convert_fq1 = f'{self.fqs_dir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.raw_convert_fq2 = f'{self.fqs_dir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        
        # rename sample
        self.refine_sample = self.sample + "_refine"
        
        # out
        self.refine_fqs_dir = os.path.abspath(f"{self.outdir}/refine_convert")
        utils.check_mkdir(self.refine_fqs_dir)
        
        self.out_fq1_file = f'{self.refine_fqs_dir}/{self.refine_sample}_S1_L001_R1_001.fastq.gz'
        self.out_fq2_file = f'{self.refine_fqs_dir}/{self.refine_sample}_S1_L001_R2_001.fastq.gz'
        self.cmd_line = f'{self.outdir}/{self.sample}_refine_cmd_line'
        
        # running parameters
        self.mem = args.mem
        self.thread = args.thread
        self.other_param = args.other_param
        self.ref_path = args.ref_path
        self.soft_path = args.soft_path

    
    @utils.add_log
    def refine_vdj(self):
        """return two sets：
        1. successfully assembled barcode sets.
        2. refined contig_id which including paired productive contig and all non-productive contig_id.
        """
        df_productive = self.raw_all_contig_anno[self.raw_all_contig_anno["productive"] == True]
        df_refine = df_productive[df_productive["is_cell"] == False]
        df_non_productive = df_refine[df_refine["productive"] == False]
        
        df_refine = df_refine.groupby("barcode").filter(lambda x: (len(x) > 4))
        df_refine = Filter_noise(df_refine, self.coeff, self.seqtype)()
        refined_contigs = set(df_refine.contig_id)
        refined_contigs |= set(df_non_productive.contig_id)

        df_productive = self.raw_filter_contig_anno[self.raw_filter_contig_anno["productive"] == True]
        raw_assembled_cells = set(df_productive.barcode)
        
        return raw_assembled_cells, refined_contigs


    @utils.add_log
    def gen_refine_fqs(self, raw_assembled_cells, refined_contigs):
        # some reads may not map to assembled contigs but used for assembling
        refined_contigs.add(None)
        read_name_set = set()

        with pysam.AlignmentFile(self.all_bam) as f:
            for read in f:
                # 记录所有被判定为细胞的contig的reads，确保不影响已经成功组装的细胞数
                if read.get_tag("CB") in raw_assembled_cells:
                    read_name_set.add(read.query_name)
                # 记录refined contig的reads
                elif read.reference_name in refined_contigs:
                    read_name_set.add(read.query_name)

        out_fq1 = xopen(self.out_fq1_file, 'w')
        out_fq2 = xopen(self.out_fq2_file, 'w')

        with pysam.FastxFile(self.raw_convert_fq1, persist=False) as fq1, \
                pysam.FastxFile(self.raw_convert_fq2, persist=False) as fq2:
            for entry1, entry2 in zip(fq1, fq2):
                header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                if header1 in read_name_set:
                    out_fq1.write(f"@{header1}\n{seq1}\n+\n{qual1}\n")
                    out_fq2.write(f'@{header2}\n{seq2}\n+\n{qual2}\n')

        out_fq1.close()
        out_fq2.close() 

    
    @utils.add_log
    def refine_assemble(self):
        """Cellranger vdj"""
        cmd = (
            f'{self.soft_path} vdj '
            f'--id={self.refine_sample} '
            f'--reference={self.ref_path} '
            f'--fastqs={self.refine_fqs_dir} '
            f'--sample={self.refine_sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )

        if self.other_param:
            cmd += (" " + self.other_param)

        with open(self.cmd_line, 'w') as f:
            f.write(cmd)
        
        cwd = os.getcwd()
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)
        # change dir back to avoid can not find '03.assemble/stat.txt' error
        os.chdir(cwd)

    def run(self):
        raw_assembled_cells, refined_contigs = self.refine_vdj()
        self.gen_refine_fqs(raw_assembled_cells, refined_contigs)
        self.refine_assemble()