import glob
import os
import pysam
import pandas as pd
import abc
import subprocess
import celescope
from Bio.Seq import Seq
from xopen import xopen

from celescope.tools.step import Step
from celescope.tools import utils


TOOLS_DIR = os.path.dirname(celescope.tools.__file__)


class NotFoundPath(Exception):
    pass


class VDJ_Mixin(Step):
    """Mixin class for cellranger"""
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title)

        self.species = args.species
        self.soft = args.soft
        self.other_param = args.other_param
        self.not_split_R2 = args.not_split_R2
        self.mem = args.mem

        if args.ref_path and args.soft_path:
            self.ref_path = args.ref_path
            self.soft_path = args.soft_path
        else:
            self.ref_path, self.soft_path = self.get_ref_path(self.soft, self.species), self.get_soft_path(self.soft)

        self.BARCODES_10X_FILE = os.path.dirname(self.soft_path) + "/lib/python/cellranger/barcodes/737K-august-2016.txt"
        self.out_fq1_file = f'{self.outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.out_fq2_file = f'{self.outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        self.barcode_correspondence_file = f'{self.outdir}/barcode_correspond.txt'

    @staticmethod
    def get_ref_path(soft, species):
        """Reutrn cellranger reference path

        :param soft: soft version
        :param species: species info
        :raises NotFoundPath
        :return ref_path
        """
        try:
            ref_path = glob.glob(f'/SGRNJ/Database/script/soft/cellranger/vdj_ref/{soft}/{species}/refdata*')
        except IndexError as e:
            raise NotFoundPath from e
        
        ref_path = [i for i in ref_path if os.path.isdir(i)][0]
        return ref_path

    @staticmethod
    def get_soft_path(soft):
        """Return cellranger soft path

        :param soft: soft version
        :raises NotFoundPath
        :return soft_path
        """
        try:
            soft_path = glob.glob(f'/SGRNJ/Database/script/soft/cellranger/cellranger-{soft}/cellranger')[0]
        except IndexError as e:
            raise NotFoundPath from e
        
        return soft_path
    
    @staticmethod
    def reversed_compl(seq):
        """Reverse complement sequence"""
        return str(Seq(seq).reverse_complement())

    @staticmethod
    def fastq_line(name, seq, qual):
        """Format fastq file

        :param name: sequence name
        :param seq: sequence
        :param qual: sequence quality
        :return: Formatted fastq file
        """
        return f'@{name}\n{seq}\n+\n{qual}\n'

    @staticmethod
    def convert_seq(sgr_barcode, sgr_umi, barcode_dict, barcodes_10X, seq2, qual2):
        """Convert each seq info of sgr format to 10X format

        :param sgr_barcode: sgr barcode info
        :param sgr_umi: sgr umi info
        :param barcode_dict: key:SGR barcode; value:10X barcode
        :param barcodes_10X: 10X barcode info
        :param seq2: r2 sequence
        :param qual2: r2 sequence quality
        :return: sequence related info in 10X format
        """
        UMI_10X_LEN = 10
        TSO = "TTTCTTATATGGG"
        if sgr_barcode in barcode_dict:
            barcode_10X = barcode_dict[sgr_barcode]
        else:
            barcode_10X = barcodes_10X.readline().strip()
            barcode_dict[sgr_barcode] = barcode_10X

        if len(sgr_umi) > UMI_10X_LEN:
            umi_10X = sgr_umi[0:UMI_10X_LEN]
        elif len(sgr_umi) < UMI_10X_LEN:
            umi_10X = sgr_umi + 'C' * (UMI_10X_LEN - len(sgr_umi))
        else:
            umi_10X = sgr_umi

        seq2_insert = 90
        seq2_cut = 60

        new_seq2_1 = seq2[0:seq2_insert]
        new_seq2_2 = seq2[seq2_cut:] 

        new_seq1 = barcode_10X + umi_10X + TSO
        new_qual1 = 'J' * len(new_seq1)
        new_qual2_1 = qual2[0:len(new_seq2_1)]
        new_qual2_2 = qual2[seq2_cut:]

        return new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2

    @utils.add_log
    def run_convert(self):
        """Convert fastq file to new file in 10X format"""
        
        barcode_dict = {}
        barcodes_10X = open(self.BARCODES_10X_FILE, 'r')
        fq_file = pysam.FastxFile(self.fq2)
        
        out_fq1 = xopen(self.out_fq1_file, 'w')
        out_fq2 = xopen(self.out_fq2_file, 'w')
        
        for entry in fq_file:
            name = entry.name
            attrs = name.split('_')
            barcode, umi = attrs[0], attrs[1]
            seq, qual = entry.sequence, entry.quality
            new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2 = self.convert_seq(barcode, umi, barcode_dict, barcodes_10X, seq, qual)

            if not self.not_split_R2:
                out_fq1.write(self.fastq_line(f'{name}_1', new_seq1, new_qual1))
                out_fq1.write(self.fastq_line(f'{name}_2', new_seq1, new_qual1))
                out_fq2.write(self.fastq_line(f'{name}_1', new_seq2_1, new_qual2_1))
                out_fq2.write(self.fastq_line(f'{name}_2', new_seq2_2, new_qual2_2))
            else:
                out_fq1.write(self.fastq_line(name, new_seq1, new_qual1))
                out_fq2.write(self.fastq_line(name, seq, qual))

        out_fq1.close()
        out_fq2.close()
        barcodes_10X.close()

        barcode_record = pd.DataFrame()
        barcode_record['sgr'] = list(barcode_dict.keys())
        barcode_record['10X'] = [barcode_dict[i] for i in barcode_dict]
        barcode_record.to_csv(self.barcode_correspondence_file, sep='\t', index=False)

    @utils.add_log
    def run_assemble(self):
        """Cellranger vdj"""
        cmd = (
            f'{self.soft_path} vdj '
            f'--id={self.sample} '
            f'--reference={self.ref_path} '
            f'--fastqs={self.fqs_dir} '
            f'--sample={self.sample} '
            f'--localcores={self.thread} '
            f'--localmem={self.mem} '
        )

        if self.other_param:
            cmd += (self.other_param)

        VDJ_Mixin.run_assemble.logger.info(cmd)
        with open(f'{self.outdir}/{self.sample}_vdj_10X.sh', 'w') as f:
            f.write(cmd)
        os.chdir(self.outdir)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run_mapping(self):
        cmd = (
            f'Rscript {TOOLS_DIR}/VDJmapping.R '
            f'--rds {self.rds} '
            f'--VDJ {self.contig} '
            f'--sample {self.sample} '
            f'--outdir {self.outdir} '
            f'--assign_file {self.assign}'
        )
        subprocess.check_call(cmd, shell=True)

    @abc.abstractmethod
    def run(self):
        pass


def get_opts_VDJ_Mixin(parser):
    parser.add_argument('--species', help='species',
                        choices=['hs', 'mmu'], default='hs')
    parser.add_argument('--soft', help='cellranger version', choices=['3.0.2', '3.1.0', '4.0.0', '6.0.0'],
                        default='4.0.0')
    parser.add_argument('--ref_path', help='reference path for cellranger')
    parser.add_argument('--soft_path', help='soft path for cellranger')
    parser.add_argument('--other_param', help='Other cellranger parameters.', default="")
    parser.add_argument('--not_split_R2', help='whether split r2',action='store_true')
    parser.add_argument('--mem', help='memory (G)', default=10)

    return parser

