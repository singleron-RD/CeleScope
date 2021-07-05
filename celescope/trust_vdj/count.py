import pysam
from celescope.tools.step import Step, s_common
from collections import defaultdict
import subprocess
import pandas as pd
from celescope.tools import utils


class Count(Step):
    """
    Features
    
    - Count umi.
    - Generate fastq file for assemble.
    
    Output
    - `count.txt` Umi types count file, sorted by umi numbers.
    - `{sample}_mapped_R1.fq` R1 reads contained barcode and umi.
    - `{sample}_mapped_R2.fq` R2 reads contained sequence.
    - `{sample}.toassemble_bc.txt` Containing barcodes that mapping to any V(D)J genes.
    """
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.sample = args.sample
        self.bam_file = args.bam_file
        self.cells = args.cells

    @utils.add_log
    def count_bam(self):
        rawbam = pysam.AlignmentFile(self.bam_file, "rb")
        header = rawbam.header
        new_bam = pysam.AlignmentFile(
            self.bam_file + ".temp", "wb", header=header)
        dic = defaultdict(set)
        for read in rawbam:
            attr = read.query_name.split('_')
            barcode = attr[0]
            umi = attr[1]
            read.set_tag(tag='CB', value=barcode, value_type='Z')
            read.set_tag(tag='UB', value=umi, value_type='Z')
            dic[barcode].add(umi)
            new_bam.write(read)
        new_bam.close()
        cmd = f'mv {self.bam_file}.temp {self.bam_file}'
        subprocess.check_call(cmd, shell=True)

        df = pd.DataFrame()
        df['barcode'] = list(dic.keys())
        df['UMI'] = [len(dic[i]) for i in list(dic.keys())]
        df = df.sort_values(by='UMI', ascending=False)
        df.to_csv(f'{self.outdir}/count.txt', sep='\t', index=False)
        return df

    @utils.add_log
    def cut_off(self):
        df = self.count_bam()

        UMI_num = int(self.cells)
        rank = UMI_num / 100
        rank_UMI = df.loc[rank, 'UMI']
        UMI_min = int(rank_UMI / 10)

        df_umi_filtered = df[df.UMI >= UMI_min]

        barcodes = df_umi_filtered.barcode.tolist()

        return barcodes


    @utils.add_log
    def filter_bam(self):
        barcodes = self.cut_off()
        bam = pysam.AlignmentFile(self.bam_file, "rb")
        new_fq1 = open(f'{self.outdir}/{self.sample}_mapped_R1.fq', 'w')
        new_fq2 = open(f'{self.outdir}/{self.sample}_mapped_R2.fq', 'w')
        R1 = set()
        R2 = set()
        for read in bam:
            attr = read.query_name.split('_')
            barcode = attr[0]
            umi = attr[1]
            quality = 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
            if barcode in barcodes:
                R1.add(f'@{read.query_name}\n{barcode}{umi}\n+\n{quality}\n')
                R2.add(f'@{read.query_name}\n{read.query_sequence}\n+\n{read.qual}\n')
            else:
                continue
        R1 = sorted(list(R1))
        R2 = sorted(list(R2))
        for r1 in R1:
            new_fq1.write(r1)
        for r2 in R2:
            new_fq2.write(r2)

        new_fq2.close()
        new_fq1.close()

        with open(f'{self.outdir}/{self.sample}.toassemble_bc.txt', 'w') as fh:
            for barcode in barcodes:
                fh.write(str(barcode) + '\n')

        Count.filter_bam.logger.info(f'get {len(barcodes)} cells for assemble')

    @utils.add_log
    def run(self):
        self.filter_bam()

@utils.add_log
def count(args):
    step_name = 'count'
    count_obj = Count(args, step_name)
    count_obj.run()


def get_opts_count(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--bam_file', help='BAM file form STAR step', required=True)
    parser.add_argument('--cells', help='expected cell number', default=3000)



            

