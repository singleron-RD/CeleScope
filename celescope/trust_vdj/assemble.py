import glob
import os
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.trust_vdj.__init__ import CHAIN, CONDA_PATH, INDEX, TOOLS_DIR


@utils.add_log
def mapping(thread, species, index_prefix, outdir, sample, fq1, fq2, barcodeRange, umiRange):
    cmd = (
        f'{CONDA_PATH}/bin/fastq-extractor -t {thread} '
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
    Assemble.process.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def trust_assemble(thread, species, outdir, sample, trimLevel=1):

    cmd = (
        f'{CONDA_PATH}/bin/trust4 -t {thread} '
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
def get_full_len_assembly(filedir, sample):
    cmd = (
        f'perl {TOOLS_DIR}/GetFullLengthAssembly.pl '
        f'{filedir}/{sample}_annot.fa > '
        f'{filedir}/{sample}_full_len.fa '
    )
    get_full_len_assembly.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def annotate(sample, thread, outdir, species):

    cmd = (
        f'{CONDA_PATH}/bin/annotator -f {INDEX}/{species}/IMGT+C.fa '
        f'-a {outdir}/{sample}_final.out '
        f'-t {thread} '
        f'-o {outdir}/{sample} '
        f'--barcode --UMI --noImpute '
        f'--readAssignment {outdir}/{sample}_assign.out '
        f'-r {outdir}/{sample}_assembled_reads.fa > {outdir}/{sample}_annot.fa'
    )
    annotate.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def fa_to_csv(outdir, sample):
    # file name
    full_len_fa = f'{outdir}/{sample}_full_len.fa'
    assign_file = f'{outdir}/{sample}_assign.out'
    # reads assignment 
    assignment = pd.read_csv(assign_file, sep='\t', header=None)
    assignment['read_barcode'] = assignment[0].apply(lambda x: x.split('_')[0])
    assignment['contig_barcode'] = assignment[1].apply(lambda x: x.split('_')[0])
    assignment['match_barcode'] = assignment[['read_barcode', 'contig_barcode']].apply(lambda x: x['read_barcode']==x['contig_barcode'], axis=1)
    assignment = assignment[assignment['match_barcode']==True]
    assignment['umi'] = assignment[0].apply(lambda x: x.split('_')[1])
    # write contig csv
    contigs = open(f'{outdir}/{sample}_contig.csv', 'w')
    # contigs.write('barcode\tis_cell\tcontig_id\thigh_confidence\tlength\tchain\tv_gene\td_gene\tj_gene\tc_gene\tfull_length\tproductive\tcdr3\tcdr3_nt\treads\tumis\traw_clonotype_id\traw_consensus_id\n')
    process_read = 0
    with pysam.FastxFile(full_len_fa) as fa:
        for read in fa:
            name = read.name
            comment = read.comment
            attrs = comment.split(' ')
            barcode = name.split('_')[0]
            is_cell = 'True'
            high_confidence = 'True'
            length = attrs[0]
            chain = attrs[2][:3]
            full_length = 'True'
            v_gene = attrs[2].split('(')[0]
            d_gene = attrs[3]
            j_gene = attrs[4].split('(')[0]
            c_gene = attrs[5]
            cdr3 = attrs[8].split('=')[1]
            cdr3_aa = 'None'
            productive = 'False'
            temp = assignment[assignment[1]==name]
            reads = str(len(temp[0].tolist()))
            umis = str(len(set(temp['umi'].tolist())))
            raw_consensus_id = 'None'
            raw_clonotype_id = 'None'
                
            string = '\t'.join([barcode, is_cell, name, high_confidence, length, chain, v_gene, d_gene, j_gene, c_gene, full_length, productive, cdr3_aa, cdr3, reads, umis, raw_clonotype_id, raw_consensus_id])
            contigs.write(f'{string}\n')
            process_read+=1
            if process_read % 10000 == 0:
                fa_to_csv.logger.info(f'Processed {process_read} contigs')

    contigs.close()

    # df = pd.read_csv(f'{outdir}/{sample}_contig.csv', sep='\t')
    # df['d_gene'] = df['d_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
    # df['c_gene'] = df['c_gene'].apply(lambda x: x.split('(')[0] if not x == '*' else 'None')
    # df['cdr3'] = df['cdr3_nt'].apply(lambda x: 'None' if "*" in str(Seq(x).translate()) or not len(x)%3==0 else str(Seq(x).translate()))
    # df['productive'] = df['cdr3'].apply(lambda x: True if not x=='None' else False)

    # df.to_csv(f'{outdir}/{sample}_contig.csv', sep=',')

    # return df

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
        self.expect_cells = int(args.expect_cells)

        # summarys
        self.match_summary = []
        self.mapping_summary = []

        # common variables
        self.chains = CHAIN[self.seqtype]

        # input
        self.match_barcodes, cell_num = utils.read_barcode_file(self.match_dir)
        del cell_num

        # output
        # dir 
        self.match_out = f'{self.outdir}/match'
        self.assemble_out = f'{self.outdir}/assemble'
        self.temp_dir = f'{self.assemble_out}/temp'
        # check dir
        utils.check_mkdir(self.match_out)
        utils.check_mkdir(self.assemble_out)
        utils.check_mkdir(self.temp_dir)

        # file
        self.match_fq1 = f'{self.match_out}/{self.sample}_matched_R1.fq'
        self.match_fq2 = f'{self.match_out}/{self.sample}_matched_R2.fq'

        self.matched_reads = 0

    @utils.add_log
    def get_match_fastq(self):
        out_fq1 = open(self.match_fq1, 'w')
        out_fq2 = open(self.match_fq2, 'w')
        
        read_dict = defaultdict(list)
        with pysam.FastxFile(self.fq2) as fq:
            for read in fq:
                attr = read.name.split('_')
                cb = attr[0]
                read_dict[cb].append(read)
            rna_cbs = [str(Seq(cb).reverse_complement()) for cb in self.match_barcodes]
            matched_cbs = set(list(read_dict.keys())).intersection(set(rna_cbs))
            if len(matched_cbs)==0:
                raise Exception('No matched barcodes found! Please check your match dir!')
            for cb in matched_cbs:
                for read in read_dict[cb]:
                    umi = attr[1]
                    qual = 'F' * len(cb + umi)
                    seq1 = f'@{read.name}\n{cb}{umi}\n+\n{qual}\n'
                    out_fq1.write(seq1)
                    out_fq2.write(str(read)+'\n')
                    matched_cbs.add(cb)
                    self.matched_reads += 1

            out_fq1.close()
            out_fq2.close()
        
        self.match_summary.append({
            'item': 'Matched Barcodes with scRNA-seq',
            'count': len(matched_cbs),
            'total_count': np.nan
        })
        self.match_summary.append({
            'item': 'Matched Reads with scRNA-seq',
            'count': self.matched_reads, 
            'total_count': np.nan
        })

        del read_dict
        
    @utils.add_log
    def process(self):
        # process all vdj
        cb_range = self.barcodeRange.split(' ')
        umi_range = self.umiRange.split(' ')

        # process single chain
        map_res = []
        
        map_index_prefix = ['bcrtcr'] + self.chains
        samples = [self.sample] * len(map_index_prefix)
        map_threads = [self.thread] * len(map_index_prefix)
        map_species = [self.species] * len(map_index_prefix)
        map_outdirs = [self.temp_dir] * len(map_index_prefix)
        map_fq1 = [self.match_fq1] * len(map_index_prefix)
        map_fq2 = [self.match_fq2] * len(map_index_prefix)
        map_cb_range = [cb_range] * len(map_index_prefix)
        map_umi_range = [umi_range] * len(map_index_prefix)
        with ProcessPoolExecutor(len(map_index_prefix)) as pool:
            for res in pool.map(mapping, map_threads, map_species, map_index_prefix, map_outdirs, samples, map_fq1, map_fq2, map_cb_range, map_umi_range):
                map_res.append(res) 

        with pysam.FastxFile(f'{self.temp_dir}/{self.sample}_bcrtcr.fq') as fl:
            self.mapping_summary.append({
                'item': 'Reads Mapped to Any V(D)J genes', 
                'count': len(list(fl)),
                'total_count': self.matched_reads
            })
        del fl

        for c in self.chains:
            with pysam.FastxFile(f'{self.temp_dir}/{self.sample}_{c}.fq') as f:
                self.mapping_summary.append({
                    'item': f'Reads Mapped to {c}', 
                    'count': len(list(f)), 
                    'total_count': self.matched_reads
                })
            del f

        # cutoff by umi num
        fl = pysam.FastxFile(f'{self.temp_dir}/{self.sample}_bcrtcr.fq')
        read_count_dict = defaultdict(int)
        umi_count_dict = defaultdict(set)
        read_dict = defaultdict(list)
        for read in fl:
            attrs = read.name.split('_')
            cb = attrs[0]
            umi = attrs[1]
            read_count_dict[cb]+=1
            umi_count_dict[cb].add(umi)
            read_dict[cb].append(read)
        barcode_list = list(read_count_dict.keys())
        read_count_list = [read_count_dict[i] for i in barcode_list]
        umi_count_list = [len(umi_count_dict[i]) for i in barcode_list]
        df = pd.DataFrame({'barcode': barcode_list, 
                            'read_count': read_count_list, 
                            'UMI': umi_count_list})
        df = df.sort_values(by='UMI', ascending=False)
        RANK = int(self.expect_cells / 100)
        rank_UMI = df.iloc[RANK, :]["UMI"]
        UMI_min = int(rank_UMI / 10)
        df["mark"] = df["UMI"].apply(
            lambda x: "CB" if (x >= UMI_min) else "UB")

        del barcode_list
        del read_count_list
        del umi_count_list
        
        df.to_csv(f'{self.assemble_out}/count.txt', sep='\t', index=False)

        df_to_split = df[df['mark']=='CB']
        df_to_split = df_to_split.sort_values(by='UMI', ascending=False)
        cell_cbs = df_to_split['barcode'].tolist()
        umi_count_l = df_to_split['UMI'].tolist()
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
            fq_ = open(f'{self.temp_dir}/temp_{i}.fq', 'w')
            cb_fa = open(f'{self.temp_dir}/temp_{i}_bc.fa', 'w')
            umi_fa = open(f'{self.temp_dir}/temp_{i}_umi.fa', 'w')
            for c in temp_cbs:
                for read in read_dict[c]:
                    name = read.name
                    umi = name.split('_')[1]
                    fq_.write(str(read)+'\n')
                    cb_fa.write(f'>{name}\n{c}\n')
                    umi_fa.write(f'>{name}\n{umi}\n')
            fq_.close()
            cb_fa.close()
            umi_fa.close()

        threads = [self.thread] * len(idx)
        temp_dirs = [self.temp_dir] * len(idx)
        temp_species = [self.species] * len(idx)
        temp_samples = [f'temp_{i}' for i in range(len(idx)-1)]
        ass = []
        with ProcessPoolExecutor(len(idx)-1) as pool:
            for res in pool.map(trust_assemble, threads, temp_species, temp_dirs, temp_samples):
                ass.append(res)
        annot = []
        with ProcessPoolExecutor(len(idx)-1) as pool:
            for res in pool.map(annotate, temp_samples, threads, temp_dirs, temp_species):
                annot.append(res) 
        gfl = []
        with ProcessPoolExecutor(len(idx)-1) as pool:
            for res in pool.map(get_full_len_assembly, temp_dirs, temp_samples):
                gfl.append(res)
        mk_csv = []
        with ProcessPoolExecutor(len(idx)-1) as pool:
            for res in pool.map(fa_to_csv, temp_dirs, temp_samples):
                mk_csv.append(res)

        # out put assemble results
        temp_outs_fa = glob.glob(f'{self.temp_dir}/temp_*_annot.fa')
        string = ' '.join(temp_outs_fa)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_annot.fa'
        Assemble.process.logger.info(cmd)
        os.system(cmd)

        temp_outs_cdr3 = glob.glob(f'{self.temp_dir}/temp_*_cdr3.out')
        string = ' '.join(temp_outs_cdr3)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_cdr3.out'
        Assemble.process.logger.info(cmd)
        os.system(cmd)

        temp_out_assembled_reads = glob.glob(f'{self.temp_dir}/temp_*_assembled_reads.fa')
        string = ' '.join(temp_out_assembled_reads)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_assembled_reads.fa'
        Assemble.process.logger.info(cmd)
        os.system(cmd)

        temp_out_reads_assign = glob.glob(f'{self.temp_dir}/temp_*_assign.out')
        string = ' '.join(temp_out_reads_assign)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_assign.out'
        Assemble.process.logger.info(cmd)
        os.system(cmd) 

        temp_contigs = glob.glob(f'{self.temp_dir}/temp_*_contig.csv')
        string = ' '.join(temp_contigs)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_contig.csv'
        Assemble.process.logger.info(cmd)
        os.system(cmd)

        # get_full_len_assembly(self.assemble_out, self.sample)      

        temp_out_full_len = glob.glob(f'{self.temp_dir}/temp_*_full_len.fa')
        string = ' '.join(temp_out_full_len)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_full_len.fa'
        Assemble.process.logger.info(cmd)
        os.system(cmd) 

        stat_file = self.outdir + '/stat.txt'
        all_summary = self.match_summary + self.mapping_summary
        sum_df = pd.DataFrame(all_summary, columns=['item', 'count', 'total_count'])
        utils.gen_stat(sum_df, stat_file)   

    def run(self):
        self.get_match_fastq()
        self.process()
        os.system(f'rm -rf {self.temp_dir}')
        self.clean_up()

@utils.add_log
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


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
    parser.add_argument('--expect_cells', help='Expected Cells number, INT.', default=3000)
