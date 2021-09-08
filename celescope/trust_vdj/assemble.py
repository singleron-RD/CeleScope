import os
from re import T
import subprocess

import pysam
import numpy as np
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.trust_vdj.__init__ import INDEX, TOOLS_DIR, CHAIN, CONDA_PATH


class Assemble(Step):
    """
    Features

    - Assemble TCR/BCR seq data.

    Output

    - `03.assemble/{sample}_toassemble.fq` Reads to assemble.
    - `03.assemble/{sample}_toassemble_bc.fa` Barcodes to assemble.
    - `03.assemble/{sample}_cdr3.out` All assembled CDR3 output.
    - `03.assemble/{sample}_barcode_report.tsv` Record chain information in each barcode.
    - `03.assemble/{sample}_annot.fa` Assembled annotated contig sequences.
    - `03.assemble/{sample}_assembled_reads.fa` Assembled raw reads.
    - `03.assemble/{sample}_report.tsv` Record assembled CDR3 types and count.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species
        self.seqtype = args.seqtype
        self.barcodeRange = args.barcodeRange
        self.umiRange = args.umiRange
        self.match_dir = args.match_dir
        self.trimLevel = args.trimLevel

        # summarys
        self.cell_summary = []
        self.match_summary = []

        # common variables
        self.chains = CHAIN[self.seqtype]

        # input
        self.match_barcodes, cell_num = utils.read_barcode_file(self.match_dir)

        # output
        # dir 
        self.match_out = f'{self.outdir}/match'
        self.assemble_out = f'{self.outdir}/assemble'
        self.summarize_out = f'{self.outdir}/summarize'
        utils.check_mkdir(self.match_out)
        utils.check_mkdir(self.assemble_out)
        utils.check_mkdir(self.summarize_out)

        # file
        self.match_fq1 = f'{self.match_out}/{self.sample}_matched_R1.fq'
        self.match_fq2 = f'{self.match_out}/{self.sample}_matched_R2.fq'

    @utils.add_log
    def get_match_fastq(self):
        out_fq1 = open(self.match_fq1, 'w')
        out_fq2 = open(self.match_fq2, 'w')

        matched_cbs = set()
        self.matched_reads = 0
        
        with pysam.FastxFile(self.fq2) as fq:
            for read in fq:
                attr = read.name.split('_')
                cb = attr[0]
                umi = attr[1]
                qual = 'F' * len(cb + umi)
                reversed_cb = str(Seq(cb).reverse_complement())
                if reversed_cb in self.match_barcodes:
                    seq1 = f'@{read.name}\n{cb}{umi}\n+\n{qual}\n'
                    out_fq1.write(seq1)
                    out_fq2.write(str(read)+'\n')
                    matched_cbs.add(reversed_cb)
                    self.matched_reads += 1

            out_fq1.close()
            out_fq2.close()
        
        self.match_summary.append({
            'item': 'Matched Barcodes',
            'count': len(matched_cbs),
            'total_count': np.nan
        })
        self.match_summary.append({
            'item': 'Matched Reads',
            'count': self.matched_reads, 
            'total_count': np.nan
        })

    @utils.add_log
    def mapping(self):
        # mapping all vdj
        cb_range = self.barcodeRange.split(' ')
        umi_range = self.umiRange.split(' ')
        cmd = (
            f'{CONDA_PATH}/bin/fastq-extractor -t {self.thread} '
            f'-f {INDEX}/{self.species}/bcrtcr.fa '
            f'-o {self.assemble_out}/{self.sample} '
            f'--barcodeStart {int(cb_range[0])} '
            f'--barcodeEnd {int(cb_range[1])} '
            f'--umiStart {int(umi_range[0])} '
            f'--umiEnd {int(umi_range[1])} '
            f'-u {self.match_fq2} '
            f'--barcode {self.match_fq1} '
            f'--UMI {self.match_fq1} ' 
        )
        subprocess.check_call(cmd, shell=True)

        fl = pysam.FastxFile(f'{self.assemble_out}/{self.sample}.fq')
        self.cell_summary.append({
            'item': f'Reads Mapped to Any V(D)J genes', 
            'count': len(list(fl)),
            'total_count': self.matched_reads
        })
        del fl
        # mapping single chain
        temp_dir = f'{self.assemble_out}/temp'
        utils.check_mkdir(temp_dir)
        for c in self.chains:
            cmd = (
                f'{CONDA_PATH}/bin/fastq-extractor -t {self.thread} '
                f'-f {INDEX}/{self.species}/{c}.fa '
                f'-o {temp_dir}/{c}_toassemble '
                f'--barcodeStart {int(cb_range[0])} '
                f'--barcodeEnd {int(cb_range[1])} '
                f'--umiStart {int(umi_range[0])} '
                f'--umiEnd {int(umi_range[1])} '
                f'-u {self.match_fq2} '
                f'--barcode {self.match_fq1} '
                f'--UMI {self.match_fq1} '
            )
            Assemble.mapping.logger.info(cmd)
            subprocess.check_call(cmd, shell=True)

            f = pysam.FastxFile(f'{temp_dir}/{c}_toassemble.fq')
            self.cell_summary.append({
                'item': f'Reads Mapped to {c}', 
                'count': len(list(f)), 
                'total_count': self.matched_reads
            })
            del f
        os.system('rm -rf {temp_dir}')

    @utils.add_log
    def run_assemble(self):

        cmd = (
            f'{CONDA_PATH}/bin/trust4  -t {self.thread} '
            f'-f {INDEX}/{self.species}/{self.seqtype}.fa '
            f'-o {self.assemble_out}/{self.sample} '
            f'-u {self.assemble_out}/{self.sample}.fq '
            f'--barcode {self.assemble_out}/{self.sample}_bc.fa '
            f'--UMI {self.assemble_out}/{self.sample}_umi.fa '
            f'--trimLevel {self.trimLevel}'
        )
        Assemble.run_assemble.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def annotate(self):

        cmd = (
            f'{CONDA_PATH}/bin/annotator -f {INDEX}/{self.species}/IMGT+C.fa '
            f'-a {self.assemble_out}/{self.sample}_final.out '
            f'-t {self.sample} '
            f'-o {self.assemble_out}/{self.sample} '
            f'--barcode --UMI '
            f'--readAssignment {self.assemble_out}/{self.sample}_assign.out '
            f'-r {self.assemble_out}/{self.sample}_assembled_reads.fa > {self.assemble_out}/{self.sample}_annot.fa'
        )
        Assemble.annotate.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def get_full_len_assembly(self):
        cmd = (
            f'perl {TOOLS_DIR}/GetFullLengthAssembly.pl '
            f'{self.assemble_out}/{self.sample}_annot.fa > '
            f'{self.assemble_out}/{self.sample}_full_len.fa '
        )
        Assemble.get_full_len_assembly.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    def run(self):
        if not os.path.exists(self.match_fq2):
            self.get_match_fastq()
        self.mapping()
        self.run_assemble()
        self.annotate()
        self.get_full_len_assembly()

@utils.add_log
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq2', help='R2 reads matched with scRNA-seq.', required=True)
        parser.add_argument('--match_dir', help='match scRNA-seq directory', required=True)

    parser.add_argument('--species', help='species', choices=["hg19", "hg38", "GRCm38"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='barcode range in fq1, INT INT CHAR', default='0 23 +') 
    parser.add_argument('--umiRange', help='umi range in fq1, INT INT CHAR', default='24 -1 +')
    parser.add_argument('--trimLevel', help='INT: 0: no trim; 1: trim low quality; 2: trim unmatched', default=1)   








