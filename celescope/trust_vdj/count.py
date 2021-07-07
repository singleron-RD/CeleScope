import pysam
from celescope.tools.step import Step, s_common
from collections import defaultdict
import subprocess
import pandas as pd
from celescope.tools import utils
from itertools import groupby

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
        
        self.match_bool = True
        if (not args.match_dir) or (args.match_dir == "None"):
            self.match_bool = False
        if self.match_bool:
            self.match_cell_barcodes, _match_cell_number = utils.read_barcode_file(
                args.match_dir)

    @utils.add_log
    def count_bam(self):
        rawbam = pysam.AlignmentFile(self.bam_file, "rb")
        header = rawbam.header
        new_bam = pysam.AlignmentFile(
            self.bam_file + ".temp", "wb", header=header)
        # write fq and bam
        dic = defaultdict(lambda: defaultdict(int))

        def keyfunc(x):
            return x.query_name
        for k,g in groupby(rawbam, keyfunc):
            name= k
            attrs = name.split('_')
            cb = attrs[0]
            umi = attrs[1]
            dic[cb][umi] += 1
            for read in list(g):
                read.set_tag(tag='CB', value=cb, value_type='Z')
                read.set_tag(tag='UB', value=umi, value_type='Z')
                new_bam.write(read)
        new_bam.close()
        cmd = f'mv {self.bam_file}.temp {self.bam_file}'
        subprocess.check_call(cmd, shell=True)
           
        # write UMI count        
        df_umi = open(f'{self.outdir}/{self.sample}_count.txt', 'w')
        df_umi.write(f'barcode\tUMI\tcount\n') 

        for cb in dic:
            for umi in dic[cb]:
                df_umi.write(f'{cb}\t{umi}\t{dic[cb][umi]}\n')
        df_umi.close()
                
        df = pd.read_csv(f'{self.outdir}/{self.sample}_count.txt', sep='\t')
        df = df.groupby(['barcode'], as_index=False).agg({'UMI': 'count'})
        df = df.sort_values(by='UMI', ascending=False)
        df.to_csv(f'{self.outdir}/{self.sample}_count.txt', sep='\t', index=False)
        
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
    def bam_to_fq(self):
        vdj_barcodes = self.cut_off()
        if self.match_bool:
            barcodes = list(set(vdj_barcodes).intersection(set(self.match_cell_barcodes)))
        
        bam = pysam.AlignmentFile(self.bam_file, "rb")
        
        qualities = 'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
        new_fq1 = open(f'{self.outdir}/{self.sample}_mapped_R1.fq', 'w')
        new_fq2 = open(f'{self.outdir}/{self.sample}_mapped_R2.fq', 'w')
        def keyfunc(x):
            return x.query_name, x.query_sequence, x.qual
        count_read = 0
        for k,g in groupby(bam, keyfunc):
            name, seq, qual = k
            attrs = name.split('_')
            cb = attrs[0]
            umi = attrs[1]
            if cb in barcodes:

                new_fq1.write(f'@{name}\n{cb}{umi}\n+\n{qualities}\n')
                new_fq2.write(f'@{name}\n{seq}\n+\n{qual}\n')
                count_read += 1
            else:
                continue
             
            if count_read % 100000 == 0:
                Count.bam_to_fq.logger.info(f'processed {count_read} reads')

        Count.bam_to_fq.logger.info(f'total processed {count_read} reads, done')
        
        new_fq1.close()
        with open(f'{self.outdir}/{self.sample}.toassemble_bc.txt', 'w') as fh:
            for barcode in barcodes:
                fh.write(str(barcode) + '\n')

        Count.bam_to_fq.logger.info(f'get {len(barcodes)} cells for assemble')

    @utils.add_log
    def run(self):
        self.bam_to_fq()

@utils.add_log
def count(args):
    step_name = 'count'
    count_obj = Count(args, step_name)
    count_obj.run()


def get_opts_count(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--bam_file', help='BAM file form STAR step', required=True)
        parser.add_argument('--match_dir', help='Match celescope scRNA-Seq directory.', default=None)
    parser.add_argument('--cells', help='expected cell number', default=3000)



            

