import os
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
from celescope.tracer_vdj.split_fastq import get_barcodes
from celescope.tools.barcode import *
import pysam
import pandas as pd


TRUST = '/SGRNJ03/randd/zhouxin/software/TRUST4/run-trust4'


def count_fq(fq1):
    bcs, umis, names = [], [], []
    count_df = pd.DataFrame()
    with pysam.FastxFile(fq1) as fq:
        for entry in fq:
            attr = entry.sequence
            cb = attr[:24]
            umi = attr[24:]
            name = entry.name
            bcs.append(cb)
            umis.append(umi)
            names.append(name)
    count_df['barcode'] = bcs
    count_df['UMI'] = umis
    count_df['seq_name'] = names
    
    return count_df

@utils.add_log
def match_barcodes(outdir, match_dir, Seqtype, fq1):
    annotated_bcs = get_barcodes(match_dir, Seqtype)
    bcs_df = pd.DataFrame(annotated_bcs, columns=['barcode'])
    count_df = count_fq(fq1)
    df = pd.merge(bcs_df, count_df, on='barcode', how='inner')
    seqnames = df['seq_name'].tolist()
    seqlist = open(f'{outdir}/seqlist.txt', 'w')
    for name in seqnames:
        seqlist.write(str(name) + '\n')

    count_df.to_csv(f'{outdir}/count.txt', sep='\t')
    df.to_csv(f'{outdir}/matched_count.txt', sep='\t')
    

class Trust_assemble(Step):
    """
    Features

    - Get fq file
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.match_dir = args.match_dir
        self.Seqtype = args.Seqtype
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species

    
    @utils.add_log
    def getFqfile(self):
        match_barcodes(self.outdir, self.match_dir, self.Seqtype, self.fq1)

        cmd1 = (
            f'seqtk subseq {self.fq1} {self.outdir}/seqlist.txt > {self.outdir}/{self.sample}_R1.fq'
        )
        os.system(cmd1)

        cmd2 = (
            f'seqtk subseq {self.fq2} {self.outdir}/seqlist.txt > {self.outdir}/{self.sample}_R2.fq'
        )
        os.system(cmd2)


    @utils.add_log
    def run(self):

        self.getFqfile()

        species = self.species

        if species =='Mmus':
            index_file = '/SGRNJ03/randd/zhouxin/software/TRUST4/mouse/GRCm38_bcrtcr.fa'
            ref = '/SGRNJ03/randd/zhouxin/software/TRUST4/mouse/mouse_IMGT+C.fa'
        elif species == 'Hsap':
            index_file = '/SGRNJ03/randd/zhouxin/software/TRUST4/hg38_bcrtcr.fa'
            ref = '/SGRNJ03/randd/zhouxin/software/TRUST4/human_IMGT+C.fa'
        cmd = (
            f'{TRUST} -t {self.thread} '
            f'-u {self.outdir}/{self.sample}_R2.fq '
            f'--barcode {self.outdir}/{self.sample}_R1.fq '
            f'--barcodeRange 0 23 + '
            f'-f {index_file} '
            f'--ref {ref} '
            f'-o {self.sample} --od {self.outdir}/TRUST4' 
        )

        os.system(cmd)

        os.remove(f'{self.outdir}/seqlist.txt')


@utils.add_log
def trust_assemble(args):
    step_name = 'trust_assemble'
    trust_assemble_obj = Trust_assemble(args, step_name)
    trust_assemble_obj.run()


def get_opts_trust_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='R1 reads from barcode step', required=True)
        parser.add_argument('--fq2', help='R2 reads from barcode step', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
    parser.add_argument('--Seqtype', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)        








