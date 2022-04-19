import glob
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import pysam
from Bio.Seq import Seq

from celescope.fl_vdj_TRUST4_split import trust_utils as tr
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.fl_vdj_TRUST4_split.__init__ import CHAIN


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
        self.trimLevel = args.trimLevel
        self.UMI_min = args.UMI_min

        self.chains = self._get_chain_type(self.seqtype)
        self.match_barcodes = self._get_match_barcode(self.match_dir)

        # outdir
        self.match_out = f'{self.outdir}/match'
        self.assemble_out = f'{self.outdir}/assemble'
        self.temp_dir = f'{self.assemble_out}/temp'
        self._check_outdir(self.match_out, self.assemble_out, self.temp_dir)
        self.matched_reads = 0
        self.match_fq1, self.match_fq2 = f'{self.match_out}/{self.sample}_matched_R1.fq', f'{self.match_out}/{self.sample}_matched_R2.fq'

    @staticmethod
    def reversed_compl(seq):
        return str(Seq(seq).reverse_complement())

    @staticmethod
    def _get_chain_type(seqtype):
        return CHAIN[seqtype]
    
    @staticmethod
    def _get_match_barcode(match_dir):
        match_barcodes, _ = utils.get_barcode_from_match_dir(match_dir) 
        return match_barcodes
    
    @staticmethod
    def _check_outdir(match_out, assemble_out, temp_dir):
        utils.check_mkdir(match_out)
        utils.check_mkdir(assemble_out)
        utils.check_mkdir(temp_dir)

    @staticmethod
    # UMI cut-off filter for candidate_reads
    def umi_cutoff(candidate_reads, UMI_min, assemble_out) :
        read_count_dict = defaultdict(int)
        umi_dict = defaultdict(set)
        read_barcode_dict = defaultdict(list)
        for read in candidate_reads:
            attrs = read.name.split('_')
            cb, umi = attrs[0], attrs[1]
            read_count_dict[cb] += 1
            umi_dict[cb].add(umi)
            read_barcode_dict[cb].append(read)
        barcode_list = list(read_count_dict.keys())

        df_count = pd.DataFrame({'barcode': list(read_count_dict.keys()), 
                            'read_count': [read_count_dict[i] for i in barcode_list], 
                            'UMI': [len(umi_dict[i]) for i in barcode_list]})

        df_count.sort_values(by='UMI', ascending=False, inplace=True)
        if UMI_min == "auto":
            RANK = 20
            rank_UMI = df_count.iloc[RANK, :]["UMI"]
            UMI_min = int(rank_UMI / 10)
        UMI_min = int(UMI_min)
        df_count["mark"] = df_count["UMI"].apply(
            lambda x: "CB" if (x >= UMI_min) else "UB")
        df_count.to_csv(f'{assemble_out}/count.txt', sep='\t', index=False)

        return df_count, read_barcode_dict

    @staticmethod
    # split candidate reads file for Multithreaded assembly
    def split_candidate_reads(df_count, read_barcode_dict, temp_dir):
        df_count_cb = df_count[df_count['mark']=='CB']
        df_count_cb.sort_values(by='UMI', ascending=False, inplace = True)
        cell_cbs = df_count_cb['barcode'].tolist()
        umi_count_l = df_count_cb['UMI'].tolist()
        sum_umi = sum(umi_count_l)
        threshold = int(sum_umi/4)
        umi_num = 0
        idx = []
        for i in range(len(umi_count_l)):
            umi_num += umi_count_l[i]
            umi_num_next = umi_num + umi_count_l[i+1]
            if umi_num<=threshold and umi_num_next>=threshold:
                idx.append(i+1)
                umi_num = 0
                if len(idx)==3:
                    break
        idx.insert(0, 0)
        idx.append(len(cell_cbs))
        for i in range(len(idx)-1):
            temp_cbs = cell_cbs[idx[i]:idx[i+1]]
            toassemble_fq = open(f'{temp_dir}/temp_{i}.fq', 'w')
            toassemble_bc = open(f'{temp_dir}/temp_{i}_bc.fa', 'w')
            toassemble_umi = open(f'{temp_dir}/temp_{i}_umi.fa', 'w')
            for tmp in temp_cbs:
                for read in read_barcode_dict[tmp]:
                    name = read.name
                    umi = name.split('_')[1]
                    toassemble_fq.write(str(read)+'\n')
                    toassemble_bc.write(f'>{name}\n{tmp}\n')
                    toassemble_umi.write(f'>{name}\n{umi}\n')
            toassemble_fq.close()
            toassemble_bc.close()
            toassemble_umi.close()
        
        return idx, len(idx)
    
    @staticmethod
    def Multi_Executor(threads, temp_dirs, temp_species, temp_samples, idx_len):
        assmble_res, annot_res, full_len_res, contig_res = [], [], [], []

        with ProcessPoolExecutor(idx_len-1) as pool:
            for res in pool.map(tr.trust_assemble, threads, temp_species, temp_dirs, temp_samples):
                assmble_res.append(res)

        with ProcessPoolExecutor(idx_len-1) as pool:
            for res in pool.map(tr.annotate, temp_samples, threads, temp_dirs, temp_species):
                annot_res.append(res) 

        with ProcessPoolExecutor(idx_len-1) as pool:
            for res in pool.map(tr.get_full_len_assembly, temp_dirs, temp_samples):
                full_len_res.append(res)

        with ProcessPoolExecutor(idx_len-1) as pool:
            for res in pool.map(tr.fa_to_csv, temp_dirs, temp_samples):
                contig_res.append(res)

    @utils.add_log
    def out_match_fastq(self):
        out_fq1 = open(self.match_fq1, 'w')
        out_fq2 = open(self.match_fq2, 'w')
        read_dict = defaultdict(list)

        with pysam.FastxFile(self.fq2) as fq:
            for read in fq:
                attr = read.name.split('_')
                cb = attr[0]
                read_dict[cb].append(read)
            rna_cbs = [self.reversed_compl(cb) for cb in self.match_barcodes]
            # matched_cbs = set(read_dict.keys())
            matched_cbs = set(read_dict.keys()).intersection(set(rna_cbs))
            
            assert len(matched_cbs) != 0
            for cb in matched_cbs:
                for read in read_dict[cb]:
                    umi = read.name.split('_')[1]
                    qual = 'F' * len(cb + umi)
                    seq1 = f'@{read.name}\n{cb}{umi}\n+\n{qual}\n'
                    out_fq1.write(seq1)
                    out_fq2.write(str(read)+'\n')
                    # matched_cbs.add(cb)
                    self.matched_reads += 1
            out_fq1.close()
            out_fq2.close()

        return matched_cbs

    @utils.add_log
    def VDJ_process(self):
        cb_range = self.barcodeRange.split(' ')
        umi_range = self.umiRange.split(' ')
        # candidate read extraction
        map_res = []
        map_index_prefix = ['bcrtcr'] + self.chains
        _map_len = len(map_index_prefix)
        samples = [self.sample] * _map_len
        map_threads = [self.thread] * _map_len
        map_species = [self.species] * _map_len
        map_outdirs = [self.temp_dir] * _map_len
        map_fq1 = [self.match_fq1] * _map_len
        map_fq2 = [self.match_fq2] * _map_len
        map_cb_range = [cb_range] * _map_len
        map_umi_range = [umi_range] * _map_len

        with ProcessPoolExecutor(_map_len) as pool:
            for res in pool.map(tr.extract_candidate_reads, map_threads, map_species, map_index_prefix, map_outdirs, samples, map_fq1, map_fq2, map_cb_range, map_umi_range):
                map_res.append(res)

        candidate_reads = pysam.FastxFile(f'{self.temp_dir}/{self.sample}_bcrtcr.fq')
        df_count, read_barcode_dict = self.umi_cutoff(candidate_reads, self.UMI_min, self.assemble_out)
        _, idx_len = self.split_candidate_reads(df_count, read_barcode_dict, self.temp_dir)

        threads = [self.thread] * idx_len
        temp_dirs = [self.temp_dir] * idx_len
        temp_species = [self.species] * idx_len
        temp_samples = [f'temp_{i}' for i in range(idx_len-1)]
        # multithread-run for assembly and annotation
        self.Multi_Executor(threads, temp_dirs, temp_species, temp_samples, idx_len)

    @utils.add_log
    def merge_file(self):
        temp_outs_fa = glob.glob(f'{self.temp_dir}/temp_*_annot.fa')
        string = ' '.join(temp_outs_fa)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_annot.fa'
        Assemble.merge_file.logger.info(cmd)
        os.system(cmd)

        temp_outs_cdr3 = glob.glob(f'{self.temp_dir}/temp_*_cdr3.out')
        string = ' '.join(temp_outs_cdr3)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_cdr3.out'
        Assemble.merge_file.logger.info(cmd)
        os.system(cmd)

        temp_out_assembled_reads = glob.glob(f'{self.temp_dir}/temp_*_assembled_reads.fa')
        string = ' '.join(temp_out_assembled_reads)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_assembled_reads.fa'
        Assemble.merge_file.logger.info(cmd)
        os.system(cmd)

        temp_out_reads_assign = glob.glob(f'{self.temp_dir}/temp_*_assign.out')
        string = ' '.join(temp_out_reads_assign)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_assign.out'
        Assemble.merge_file.logger.info(cmd)
        os.system(cmd) 

        temp_contigs = glob.glob(f'{self.temp_dir}/temp_*_contig.csv')
        string = ' '.join(temp_contigs)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_contig.csv'
        Assemble.merge_file.logger.info(cmd)
        os.system(cmd)     

        temp_out_full_len = glob.glob(f'{self.temp_dir}/temp_*_full_len.fa')
        string = ' '.join(temp_out_full_len)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_full_len.fa'
        Assemble.merge_file.logger.info(cmd)
        os.system(cmd) 

    def gen_all_contig_fasta(self):
        os.system(f'mkdir {self.outdir}/../04.summarize')
        full_len_fa = f'{self.assemble_out}/{self.sample}_full_len.fa'
        all_fa = open(f'{self.outdir}/../04.summarize/{self.sample}_all_contig.fasta','w')
        with pysam.FastxFile(full_len_fa) as fa:
            for read in fa: 
                name = read.name
                barcode = name.split('_')[0]
                sequence = read.sequence
                all_fa.write('>' + self.reversed_compl(barcode) + '_' + name.split('_')[1] + '\n' + sequence + '\n')    
        all_fa.close()
    
    def gen_report(self, matched_cbs):
        tr.get_trust_report(self.assemble_out,self.sample)
        tr.get_bc_report(self.assemble_out, self.sample)
        tr.get_bcfilter_report(self.assemble_out)

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

        with pysam.FastxFile(f'{self.temp_dir}/{self.sample}_bcrtcr.fq') as f:
            self.add_metric(
                name = 'Reads Mapped to Any V(D)J genes', 
                value = len(list(f)),
                total = self.matched_reads,
                help_info = "Fraction of reads that partially or wholly map to any germline V(D)J gene segment"
            )

        for _chain in self.chains:
            with pysam.FastxFile(f'{self.temp_dir}/{self.sample}_{_chain}.fq') as f:
                self.add_metric(
                    name = f'Reads Mapped to {_chain}', 
                    value = len(list(f)), 
                    total = self.matched_reads,
                    help_info = f"Fraction of reads that map partially or wholly to a germline {_chain} gene segment."
                )

    def run(self):
        matched_cbs = self.out_match_fastq()
        self.VDJ_process()
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

    parser.add_argument('--species', help='Species name and version.', choices=["hg19", "hg38", "GRCm38"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='Barcode range in fq1, INT INT CHAR.', default='0 23 +') 
    parser.add_argument('--umiRange', help='UMI range in fq1, INT INT CHAR.', default='24 -1 +')
    parser.add_argument('--trimLevel', help='INT: 0: no trim; 1: trim low quality; 2: trim unmatched.', default=1)
    parser.add_argument('--UMI_min', help='UMI number, INT.', default='auto')
