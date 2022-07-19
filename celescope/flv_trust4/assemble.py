import subprocess
from collections import defaultdict
from multiprocessing import Pool
import math
import random

import pandas as pd
import numpy as np
import pysam

from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.flv_trust4.__init__ import CHAIN, REF_DIR, TOOLS_DIR

# split fastq into const N_CHUNK. Different N_CHUNK will cause different results as discussed in 
# https://github.com/liulab-dfci/TRUST4/issues/75
SPLIT_N_CHUNKS = 4


class Assemble(Step):
    """
    ## Features

    - TRUST4 does not use multi-processing when assembling. By default, the candidate reads are split into 4 chunks to speed up.

    - Keep only full-length contigs.

    ## Output
    - `03.assemble/assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.

    - `03.assemble/assemble/{sample}_assign.out` Read assignment results.

    - `03.assemble/assemble/{sample}_assembled_reads.fa` Assembled raw reads sequence.

    - `03.assemble/assemble/{sample}_annotate.fa` Assembled annotated contig sequences info.

    - `03.assemble/assemble/{sample}_report.tsv` Record assembled CDR3 types, read count and proportion.

    - `03.assemble/assemble/{sample}_filter_report.tsv` Filter non-functional CDR3.

    - `03.assemble/assemble/{sample}_barcode_report.tsv` Record chain information for each barcode.
    
    - `03.assemble/assemble/{sample}_barcode_filter_report.tsv` If --diffuseFrac provided. Filter low abundance barcode.
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.ref = args.ref
        self.seqtype = args.seqtype
        self.barcodeRange = args.barcodeRange
        self.umiRange = args.umiRange
        self.candidate_fq = args.candidate_fq

        if args.not_split:
            self._n_chunk = 1
        else:
            self._n_chunk = SPLIT_N_CHUNKS
        self._chains = CHAIN[self.seqtype]
        # total thread may be higher than argument thread
        self._single_thread = math.ceil(self.thread / SPLIT_N_CHUNKS)

        # outdir
        self.assemble_outdir = f'{self.outdir}/assemble'
        self.temp_outdir = f'{self.assemble_outdir}/temp'
        for d in [self.assemble_outdir, self.temp_outdir]:
            utils.check_mkdir(dir_name=d)

        self.temp_outdir_list = [self.temp_outdir] * self._n_chunk
        self.temp_ref_list = [self.ref] * self._n_chunk
        self.temp_name_list = [f'temp_{i}' for i in range(self._n_chunk)]
        self.single_thread_list = [self._single_thread] * self._n_chunk

    def run(self):
        self.split_candidate_reads()
        self.run_assemble()
        self.run_annotate()
        self.merge_file()
        self.gen_report()

    @utils.add_log
    def split_candidate_reads(self):
        """
        split original candidate reads(_bcrtcr.fq) by barcode into N_CHUNK files
        """
        read_count_dict, umi_dict = defaultdict(int), defaultdict(set)
        with pysam.FastxFile(self.candidate_fq) as f:
            for read in f:
                attrs = read.name.split('_')
                cb, umi = attrs[0], attrs[1]
                read_count_dict[cb] += 1
                umi_dict[cb].add(umi)
        barcode_list = list(read_count_dict.keys())
        df_count = pd.DataFrame({'barcode': barcode_list, 
                            'read_count': [read_count_dict[i] for i in barcode_list], 
                            'UMI': [len(umi_dict[i]) for i in barcode_list]})
        df_count.sort_values(by='UMI', ascending=False, inplace=True)
        df_count.to_csv(f'{self.assemble_outdir}/count.txt', sep='\t', index=False)

        del read_count_dict
        del umi_dict

        # Split barcode list into equal size chunks.
        random.shuffle(barcode_list)
        split_barcodes = np.array_split(barcode_list, self._n_chunk)
        split_barcodes = [set(i) for i in split_barcodes]

        fq_list = [open(f'{self.temp_outdir}/temp_{i}.fq','w') for i in range(self._n_chunk)]
        bc_list = [open(f'{self.temp_outdir}/temp_{i}_bc.fa','w') for i in range(self._n_chunk)]
        umi_list = [open(f'{self.temp_outdir}/temp_{i}_umi.fa','w') for i in range(self._n_chunk)]

        with pysam.FastxFile(self.candidate_fq) as f:
            for read in f:
                name = read.name
                bc, umi = name.split('_')[0], name.split('_')[1]
                for i in range(self._n_chunk):
                    if bc in split_barcodes[i]:
                        fq_list[i].write(str(read) + '\n')
                        bc_list[i].write(f'>{name}\n{bc}\n')
                        umi_list[i].write(f'>{name}\n{umi}\n')
                        break

        for i in range(self._n_chunk):
            fq_list[i].close()
            bc_list[i].close()
            umi_list[i].close()
    
    @utils.add_log
    def run_assemble(self):
        """
        run assemble for each chunk
        """
        with Pool(self._n_chunk) as pool:
            pool.starmap(
                Assemble.assemble, 
                zip(self.temp_ref_list, self.temp_outdir_list, self.temp_name_list, self.single_thread_list)
            )

    @staticmethod
    @utils.add_log
    def assemble(ref, outdir, sample, single_thread, trimLevel=1):

        cmd = (
            f'trust4 -t {single_thread} '
            f'-f {REF_DIR}/{ref}/bcrtcr.fa '
            f'-o {outdir}/{sample} '
            f'-u {outdir}/{sample}.fq '
            f'--barcode {outdir}/{sample}_bc.fa '
            f'--UMI {outdir}/{sample}_umi.fa '
            f'--trimLevel {trimLevel} '
            '2>&1 '
        )
        Assemble.assemble.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def run_annotate(self):
        """
        run annotate for each chunk
        """
        with Pool(self._n_chunk) as pool:
            pool.starmap(Assemble.annotate, zip(self.temp_name_list, self.temp_outdir_list, self.temp_ref_list, self.single_thread_list))

    @staticmethod
    @utils.add_log
    def annotate(name, outdir, ref, single_thread):

        cmd = (
            f'annotator -f {REF_DIR}/{ref}/IMGT+C.fa '
            f'-a {outdir}/{name}_final.out '
            f'-t {single_thread} '
            f'-o {outdir}/{name} '
            f'--barcode --UMI --noImpute '
            f'--readAssignment {outdir}/{name}_assign.out '
            f'-r {outdir}/{name}_assembled_reads.fa > {outdir}/{name}_annotate.fa '
        )
        Assemble.annotate.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def merge_file(self):
        """
        merge files from all threads into one file.
        """

        file_suffixes = [
            'annotate.fa',
            'cdr3.out',
            'assembled_reads.fa',
            'assign.out',
        ]

        for file_suffix in file_suffixes:
            file_path = [f'{self.temp_outdir}/temp_{thread}_{file_suffix}' for thread in range(self._n_chunk)]
            file_path_str = ' '.join(file_path)
            cmd = f'cat {file_path_str} > {self.assemble_outdir}/{self.sample}_{file_suffix}'
            self.merge_file.logger.info(cmd)
            self.debug_subprocess_call(cmd)

    def gen_report(self):
        Assemble.get_trust_report(self.assemble_outdir, self.sample)
        Assemble.filter_trust_report(self.assemble_outdir, self.sample)
        Assemble.get_bc_report(self.assemble_outdir, self.sample)
        Assemble.get_bcfilter_report(self.assemble_outdir, self.sample)
        Assemble.get_airr_file(self.assemble_outdir, self.sample)


    @staticmethod
    @utils.add_log
    def get_trust_report(filedir, sample):
        cmd = (
            f'perl {TOOLS_DIR}/trust-simplerep.pl '
            f'{filedir}/{sample}_cdr3.out > '
            f'{filedir}/{sample}_report.tsv '
            '2>&1 '
        )
        Assemble.get_trust_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    @utils.add_log
    def filter_trust_report(filedir, sample):
        cmd = f''' awk '$4!~"_" && $4!~"?"' {filedir}/{sample}_report.tsv > {filedir}/{sample}_filter_report.tsv 2>&1'''
        Assemble.filter_trust_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    @utils.add_log
    def get_bc_report(filedir, sample):
        cmd = (
            f'perl {TOOLS_DIR}/trust-barcoderep.pl '
            f'{filedir}/{sample}_cdr3.out > '
            f'{filedir}/{sample}_barcode_report.tsv '
            '2>&1 ' 
        )
        Assemble.get_bc_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    @utils.add_log
    def get_bcfilter_report(filedir, sample):
        cmd = (
            f'python {TOOLS_DIR}/barcoderep-filter.py '
            f'-b {filedir}/{sample}_barcode_report.tsv > '
            f'{filedir}/{sample}_barcode_filter_report.tsv '
            '2>&1 '
        )
        Assemble.get_bcfilter_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    
    @staticmethod
    @utils.add_log
    def get_airr_file(filedir, sample):
        cmd = (
            f'perl {TOOLS_DIR}/trust-airr.pl '
            f'{filedir}/{sample}_barcode_report.tsv '
            f'{filedir}/{sample}_annotate.fa --format barcoderep > '
            f'{filedir}/airr_rearrangement.tsv '
            '2>&1 '
        )
        Assemble.get_airr_file.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)


@utils.add_log
def assemble(args):
    assemble_obj = Assemble(args)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--candidate_fq', help='Candidate fastq file from mapping step', required=True)

    parser.add_argument('--not_split', help='do not split reads into chunks', action='store_true')
    parser.add_argument('--ref', help='reference name', choices=["hg19", "hg38", "GRCm38", "other"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='Barcode range in fq1, INT INT CHAR.', default='0 23 +') 
    parser.add_argument('--umiRange', help='UMI range in fq1, INT INT CHAR.', default='24 -1 +')
