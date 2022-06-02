import subprocess
from collections import defaultdict
from multiprocessing import Pool
import math
import random

import pandas as pd
import numpy as np
import pysam
from Bio.Seq import Seq

from celescope.fl_vdj_TRUST4_split import trust_utils as tr
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.fl_vdj_TRUST4_split.__init__ import CHAIN

# split fastq into const N_CHUNK. Different N_CHUNK will cause different results as discussed in 
# https://github.com/liulab-dfci/TRUST4/issues/75
N_CHUNK = 4

class Assemble(Step):
    """
    ## Features

    - Assemble TCR/BCR seq data.

    ## Output
    - `03.assemble/match/{sample}_matched_R1.fq` New R1 reads matched with scRNA-seq
    - `03.assemble/match/{sample}_matched_R1.fq` New R2 reads matched with scRNA-seq

    - `03.assemble/assemble/{sample}_cdr3.out` All assembled CDR3 output and gene information.
    - `03.assemble/assemble/{sample}_assign.out` read assignment results.
    - `03.assemble/assemble/{sample}_assembled_reads.fa` Assembled raw reads.
    - `03.assemble/assemble/{sample}_annot.fa` Assembled annotated contig sequences.
    - `03.assemble/assemble/{sample}_full_len.fa` Assembled full length contig sequences.
    - `03.assemble/assemble/report.out` Record assembled CDR3 types, read count and proportion of read count.
    - `03.assemble/assemble/barcoderep.tsv` Record chain information for each barcode.
    - `03.assemble/assemble/barcoderepfl.tsv` Record chain information for each barcode(preliminary filter).
    """
    def __init__(self, args, display_title=None):
        Step.__init__(self, args, display_title=display_title)

        self.fq2 = args.fq2
        self.species = args.species
        self.seqtype = args.seqtype
        self.barcodeRange = args.barcodeRange
        self.umiRange = args.umiRange
        self.match_dir = args.match_dir

        self.chains = self._get_chain_type(self.seqtype)
        self.match_barcodes, _ = utils.get_barcode_from_match_dir(self.match_dir)
        # total thread may be higher than argument thread
        self.single_thread = math.ceil(self.thread / N_CHUNK)


        # outdir
        self.match_outdir = f'{self.outdir}/match'
        self.assemble_outdir = f'{self.outdir}/assemble'
        self.temp_outdir = f'{self.assemble_outdir}/temp'
        for d in [self.match_outdir, self.assemble_outdir, self.temp_outdir]:
            utils.check_mkdir(dir_name=d)
        self.matched_reads = 0

        # output files
        self.match_fq1, self.match_fq2 = f'{self.match_outdir}/{self.sample}_matched_R1.fq', \
            f'{self.match_outdir}/{self.sample}_matched_R2.fq'
        self.candidate_reads = f'{self.temp_outdir}/{self.sample}_bcrtcr.fq'

    @staticmethod
    def reversed_compl(seq):
        return str(Seq(seq).reverse_complement())

    @staticmethod
    def _get_chain_type(seqtype):
        return CHAIN[seqtype]

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
        split_barcodes = np.array_split(barcode_list, N_CHUNK)
        split_barcodes = [set(i) for i in split_barcodes]

        fq_list = [open(f'{self.temp_outdir}/temp_{i}.fq','w') for i in range(N_CHUNK)]
        bc_list = [open(f'{self.temp_outdir}/temp_{i}_bc.fa','w') for i in range(N_CHUNK)]
        umi_list = [open(f'{self.temp_outdir}/temp_{i}_umi.fa','w') for i in range(N_CHUNK)]

        with pysam.FastxFile(self.candidate_reads) as f:
            for read in f:
                name = read.name
                bc, umi = name.split('_')[0], name.split('_')[1]
                for i in range(N_CHUNK):
                    if bc in split_barcodes[i]:
                        fq_list[i].write(str(read) + '\n')
                        bc_list[i].write(f'>{name}\n{bc}\n')
                        umi_list[i].write(f'>{name}\n{umi}\n')
                        break

        for i in range(N_CHUNK):
            fq_list[i].close()
            bc_list[i].close()
            umi_list[i].close()
    
    @utils.add_log
    def assemble(self):
        """
        trust_assemble
        annotate
        get_full_len_assembly
        fa_to_csv
        """

        temp_outdirs = [self.temp_outdir] * N_CHUNK
        temp_species = [self.species] * N_CHUNK
        temp_samples = [f'temp_{i}' for i in range(N_CHUNK)]
        single_threads = [self.single_thread] * N_CHUNK

        with Pool(N_CHUNK) as pool:
            pool.starmap(tr.trust_assemble, zip(temp_species, temp_outdirs, temp_samples, single_threads))


        with Pool(N_CHUNK) as pool:
            pool.starmap(tr.annotate, zip(temp_samples, temp_outdirs, temp_species, single_threads))


        with Pool(N_CHUNK) as pool:
            pool.starmap(tr.get_full_len_assembly, zip(temp_outdirs, temp_samples))


        with Pool(N_CHUNK) as pool:
            pool.starmap(tr.fa_to_csv, zip(temp_outdirs, temp_samples))


    @utils.add_log
    def out_match_fastq(self):
        out_fq1 = open(self.match_fq1, 'w')
        out_fq2 = open(self.match_fq2, 'w')

        matched_cbs, rna_cbs = set(), {self.reversed_compl(cb) for cb in self.match_barcodes}
        with pysam.FastxFile(self.fq2) as fq:
            for read in fq:
                attr = read.name.split('_')
                cb = attr[0]
                umi = attr[1]
                qual = 'F' * len(cb + umi)
                seq1 = f'@{read.name}\n{cb}{umi}\n+\n{qual}\n'
                if cb in rna_cbs:
                    out_fq1.write(seq1)
                    out_fq2.write(str(read)+'\n')
                    self.matched_reads += 1
                    matched_cbs.add(cb)
        out_fq1.close()
        out_fq2.close()

        return matched_cbs

    @utils.add_log
    def extract_chain_reads(self):
        """
        bcrtcr.fq + {chains}.fq
        """
        cb_range = self.barcodeRange.split(' ')
        umi_range = self.umiRange.split(' ')

        map_index_prefix = ['bcrtcr'] + self.chains
        _map_len = len(map_index_prefix)
        samples = [self.sample] * _map_len
        map_species = [self.species] * _map_len
        map_outdirs = [self.temp_outdir] * _map_len
        map_fq1 = [self.match_fq1] * _map_len
        map_fq2 = [self.match_fq2] * _map_len
        map_cb_range = [cb_range] * _map_len
        map_umi_range = [umi_range] * _map_len
        map_n_thread = [self.single_thread] * _map_len

        with Pool(_map_len) as pool:
            pool.starmap(tr.extract_candidate_reads, 
            zip(map_species, map_index_prefix, map_outdirs, samples, map_fq1, map_fq2, map_cb_range, map_umi_range, map_n_thread))       

    @utils.add_log
    def merge_file(self):
        """
        merge files from all threads into one file.
        """

        file_suffixes = [
            'annot.fa',
            'cdr3.out',
            'assembled_reads.fa',
            'assign.out',
            'contig.csv',
            'full_len.fa',
        ]

        for file_suffix in file_suffixes:
            file_path = [f'{self.temp_outdir}/temp_{thread}_{file_suffix}' for thread in range(N_CHUNK)]
            file_path_str = ' '.join(file_path)
            cmd = f'cat {file_path_str} > {self.assemble_outdir}/{self.sample}_{file_suffix}'
            self.merge_file.logger.info(cmd)
            subprocess.check_call(cmd, shell=True)


    def gen_all_contig_fasta(self):
        utils.check_mkdir(f'{self.outdir}/../04.summarize')
        full_len_fa = f'{self.assemble_outdir}/{self.sample}_full_len.fa'
        all_fa = open(f'{self.outdir}/../04.summarize/{self.sample}_all_contig.fasta','w')
        with pysam.FastxFile(full_len_fa) as fa:
            for read in fa: 
                name = read.name
                barcode = name.split('_')[0]
                sequence = read.sequence
                all_fa.write('>' + self.reversed_compl(barcode) + '_' + name.split('_')[1] + '\n' + sequence + '\n')    
        all_fa.close()
    
    def gen_report(self, matched_cbs):
        tr.get_trust_report(self.assemble_outdir,self.sample)
        tr.filter_trust_report(self.assemble_outdir)
        tr.get_bc_report(self.assemble_outdir, self.sample)
        tr.get_bcfilter_report(self.assemble_outdir)

        self.add_metric(
            name="Matched Barcodes with scRNA-seq",
            value=len(matched_cbs),
            help_info="Number of barcodes those were determined to be cells in scRNA-seq"
            )

        self.add_metric(
            name = 'Matched Reads with scRNA-seq',
            value = self.matched_reads, 
            help_info="Number of reads those were from barcodes determined to be cells in scRNA-seq"
            )

        with pysam.FastxFile(f'{self.temp_outdir}/{self.sample}_bcrtcr.fq') as f:
            self.add_metric(
                name = 'Reads Mapped to Any V(D)J genes', 
                value = len(list(f)),
                total = self.matched_reads,
                help_info = "Fraction of reads that partially or wholly map to any germline V(D)J gene segment"
            )

        for _chain in self.chains:
            with pysam.FastxFile(f'{self.temp_outdir}/{self.sample}_{_chain}.fq') as f:
                self.add_metric(
                    name = f'Reads Mapped to {_chain}', 
                    value = len(list(f)), 
                    total = self.matched_reads,
                    help_info = f"Fraction of reads that map partially or wholly to a germline {_chain} gene segment."
                )

    def run(self):
        matched_cbs = self.out_match_fastq()
        self.extract_chain_reads()
        self.split_candidate_reads()
        self.assemble()
        self.merge_file()
        self.gen_all_contig_fasta()
        self.gen_report(matched_cbs)


@utils.add_log
def assemble(args):
    with Assemble(args, display_title="Match and Mapping") as runner:
        runner.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq2', help='R2 reads matched with scRNA-seq.', required=True)
        parser.add_argument('--match_dir', help='Match scRNA-seq directory.', required=True)

    parser.add_argument('--species', help='Species name and version.', choices=["hg19", "hg38", "GRCm38", "other"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='Barcode range in fq1, INT INT CHAR.', default='0 23 +') 
    parser.add_argument('--umiRange', help='UMI range in fq1, INT INT CHAR.', default='24 -1 +')
