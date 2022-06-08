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
from celescope.flv_trust4.__init__ import CHAIN

# split fastq into const N_CHUNK. Different N_CHUNK will cause different results as discussed in 
# https://github.com/liulab-dfci/TRUST4/issues/75
SPLIT_N_CHUNKS = 4

ANNOTATE_FA_SUFFIX = 'annotate.fa'
REPORT_SUFFIX = 'report.tsv'

class Assemble(Step):
    """
    ## Features

    - Assemble TCR/BCR seq data.

    ## Output
    - `03.assemble/assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.
    - `03.assemble/assemble/{sample}_assign.out` read assignment results.
    - `03.assemble/assemble/{sample}_assembled_reads.fa` Assembled raw reads.
    - `03.assemble/assemble/{sample}_annot.fa` Assembled annotated contig sequences.
    - `03.assemble/assemble/report.tsv` Record assembled CDR3 types, read count and proportion of read count.
    - `03.assemble/assemble/barcoderep.tsv` Record chain information for each barcode.
    - `03.assemble/assemble/barcoderepfl.tsv` Record chain information for each barcode(preliminary filter).
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.ref = args.ref
        self.seqtype = args.seqtype
        self.barcodeRange = args.barcodeRange
        self.umiRange = args.umiRange
        self.match_fq1 = args.match_fq1
        self.match_fq2 = args.match_fq2

        if args.not_split:
            self._n_chunk = 1
        else:
            self._n_chunk = SPLIT_N_CHUNKS
        self._chains = CHAIN[self.seqtype]
        # total thread may be higher than argument thread
        self._single_thread = math.ceil(self.thread / SPLIT_N_CHUNKS)
        self._match_reads = self.get_slot_key(
            slot='metrics',
            step_name='barcode',
            key='Valid Matched Reads',
        )

        # outdir
        self.assemble_outdir = f'{self.outdir}/assemble'
        self.temp_outdir = f'{self.assemble_outdir}/temp'
        for d in [self.assemble_outdir, self.temp_outdir]:
            utils.check_mkdir(dir_name=d)

        # output files
        self.candidate_reads = f'{self.temp_outdir}/{self.sample}_bcrtcr.fq'

        # 
        self.temp_outdir_list = [self.temp_outdir] * SPLIT_N_CHUNKS
        self.temp_ref_list = [self.ref] * SPLIT_N_CHUNKS
        self.temp_name_list = [f'temp_{i}' for i in range(SPLIT_N_CHUNKS)]
        self.single_thread_list = [self._single_thread] * SPLIT_N_CHUNKS

    def run(self):
        self.extract_chain_reads()
        self.split_candidate_reads()
        self.assemble()
        self.merge_file()
        self.gen_report()

    @utils.add_log
    def extract_chain_reads(self):
        """
        bcrtcr.fq + {chains}.fq
        """
        cb_range = self.barcodeRange.split(' ')
        umi_range = self.umiRange.split(' ')

        map_index_prefix = ['bcrtcr'] + self._chains
        n_map = len(map_index_prefix)
        samples = [self.sample] * n_map
        map_ref = [self.ref] * n_map
        map_outdirs = [self.temp_outdir] * n_map
        map_fq1 = [self.match_fq1] * n_map
        map_fq2 = [self.match_fq2] * n_map
        map_cb_range = [cb_range] * n_map
        map_umi_range = [umi_range] * n_map
        map_n_thread = [self._single_thread] * n_map

        with Pool(n_map) as pool:
            pool.starmap(Assemble.extract_candidate_reads, 
            zip(map_ref, map_index_prefix, map_outdirs, samples, map_fq1, map_fq2, map_cb_range, map_umi_range, map_n_thread))  

    @staticmethod
    @utils.add_log
    def extract_candidate_reads(species, index_prefix, outdir, sample, fq1, fq2, barcodeRange, umiRange, single_thread):
        """
        helper function for extract reads map to index_prefix
        index_prefix can be 'bcrtcr' + chains
        """
        cmd = (
            f'fastq-extractor -t {single_thread} '
            f'-f {INDEX}/{species}/{index_prefix}.fa '
            f'-o {outdir}/{sample}_{index_prefix} '
            f'--barcodeStart {barcodeRange[0]} '
            f'--barcodeEnd {barcodeRange[1]} '
            f'--umiStart {umiRange[0]} '
            f'--umiEnd {umiRange[1]} '
            f'-u {fq2} '
            f'--barcode {fq1} '
            f'--UMI {fq1} '
            )
        extract_candidate_reads.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)  

    @utils.add_log
    def split_candidate_reads(self):
        """
        split original candidate reads(_bcrtcr.fq) by barcode into N_CHUNK files
        """
        read_count_dict, umi_dict = defaultdict(int), defaultdict(set)
        with pysam.FastxFile(self.candidate_reads) as f:
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
        split_barcodes = np.array_split(barcode_list, SPLIT_N_CHUNKS)
        split_barcodes = [set(i) for i in split_barcodes]

        fq_list = [open(f'{self.temp_outdir}/temp_{i}.fq','w') for i in range(SPLIT_N_CHUNKS)]
        bc_list = [open(f'{self.temp_outdir}/temp_{i}_bc.fa','w') for i in range(SPLIT_N_CHUNKS)]
        umi_list = [open(f'{self.temp_outdir}/temp_{i}_umi.fa','w') for i in range(SPLIT_N_CHUNKS)]

        with pysam.FastxFile(self.candidate_reads) as f:
            for read in f:
                name = read.name
                bc, umi = name.split('_')[0], name.split('_')[1]
                for i in range(SPLIT_N_CHUNKS):
                    if bc in split_barcodes[i]:
                        fq_list[i].write(str(read) + '\n')
                        bc_list[i].write(f'>{name}\n{bc}\n')
                        umi_list[i].write(f'>{name}\n{umi}\n')
                        break

        for i in range(SPLIT_N_CHUNKS):
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
                Assemble.trust_assemble, 
                zip(self.temp_ref_list, self.temp_outdir_list, self.temp_name_list, self.single_thread_list)
            )

    @staticmethod
    @utils.add_log
    def assemble(species, outdir, sample, single_thread, trimLevel=1):

        cmd = (
            f'trust4 -t {single_thread} '
            f'-f {INDEX}/{species}/bcrtcr.fa '
            f'-o {outdir}/{sample} '
            f'-u {outdir}/{sample}.fq '
            f'--barcode {outdir}/{sample}_bc.fa '
            f'--UMI {outdir}/{sample}_umi.fa '
            f'--trimLevel {trimLevel}'
        )
        trust_assemble.logger.info(cmd)
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
    def annotate(name, outdir, species, single_thread):

        cmd = (
            f'annotator -f {INDEX}/{species}/IMGT+C.fa '
            f'-a {outdir}/{name}_final.out '
            f'-t {single_thread} '
            f'-o {outdir}/{name} '
            f'--barcode --UMI --noImpute '
            f'--readAssignment {outdir}/{name}_assign.out '
            f'-r {outdir}/{name}_assembled_reads.fa > {outdir}/{name}_{ANNOTATE_FA_SUFFIX}'
        )
        annotate.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def out_match_fastq(self):
        """
        Count matched barcodes and matched reads.
        """
        matched_cbs = set()
        with pysam.FastxFile(self.match_fq2) as fq:
            for read in fq:
                cb = read.name.split('_')[0]
                self.matched_reads += 1
                matched_cbs.add(cb)

        return matched_cbs


    @utils.add_log
    def merge_file(self):
        """
        merge files from all threads into one file.
        """

        file_suffixes = [
            ANNOTATE_FA_SUFFIX,
            'cdr3.out',
            'assembled_reads.fa',
            'assign.out',
        ]

        for file_suffix in file_suffixes:
            file_path = [f'{self.temp_outdir}/temp_{thread}_{file_suffix}' for thread in range(SPLIT_N_CHUNKS)]
            file_path_str = ' '.join(file_path)
            cmd = f'cat {file_path_str} > {self.assemble_outdir}/{self.sample}_{file_suffix}'
            self.merge_file.logger.info(cmd)
            subprocess.check_call(cmd, shell=True)

    def gen_report(self):
        Assemble.get_trust_report(self.assemble_outdir, self.sample)
        Assemble.filter_trust_report(self.assemble_outdir)
        Assemble.get_bc_report(self.assemble_outdir, self.sample)
        Assemble.get_bcfilter_report(self.assemble_outdir)

        with pysam.FastxFile(f'{self.temp_outdir}/{self.sample}_bcrtcr.fq') as f:
            self.add_metric(
                name = 'Reads Mapped to Any V(D)J genes', 
                value = len(list(f)),
                total = self.matched_reads,
                help_info = "Fraction of reads that partially or wholly map to any germline V(D)J gene segment"
            )

        for _chain in self._chains:
            with pysam.FastxFile(f'{self.temp_outdir}/{self.sample}_{_chain}.fq') as f:
                self.add_metric(
                    name = f'Reads Mapped to {_chain}', 
                    value = len(list(f)), 
                    total = self.matched_reads,
                    help_info = f"Fraction of reads that map partially or wholly to a germline {_chain} gene segment."
                )

    @staticmethod
    @utils.add_log
    def get_trust_report(filedir, sample):
        cmd = (
            f'perl {TOOLS_DIR}/trust-simplerep.pl '
            f'{filedir}/{sample}_cdr3.out > '
            f'{filedir}/report.tsv'
        )
        get_trust_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    @utils.add_log
    def filter_trust_report(filedir):
        cmd = f''' awk '$4!~"_" && $4!~"?"' {filedir}/report.tsv > {filedir}/reportfl.tsv '''
        filter_trust_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    @utils.add_log
    def get_bc_report(filedir, sample):
        cmd = (
            f'perl {TOOLS_DIR}/trust-barcoderep.pl '
            f'{filedir}/{sample}_cdr3.out > '
            f'{filedir}/barcoderep.tsv ' 
        )
        get_bc_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    @utils.add_log
    def get_bcfilter_report(filedir):
        cmd = (
            f'python {TOOLS_DIR}/barcoderep-filter.py '
            f'-b {filedir}/barcoderep.tsv > '
            f'{filedir}/barcoderepfl.tsv '
        )
        get_bcfilter_report.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

@utils.add_log
def assemble(args):
    with Assemble(args, display_title="Match and Mapping") as runner:
        runner.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--candidate_fq', help='Candidate fastq file from mapping step', required=True)

    parser.add_argument('--not_split', help='do not split reads into chunks', action='store_true')
    parser.add_argument('--ref', help='reference name', choices=["hg19", "hg38", "GRCm38", "other"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='Barcode range in fq1, INT INT CHAR.', default='0 23 +') 
    parser.add_argument('--umiRange', help='UMI range in fq1, INT INT CHAR.', default='24 -1 +')