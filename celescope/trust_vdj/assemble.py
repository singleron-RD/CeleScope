import os
from re import T
import subprocess
import glob
from collections import defaultdict

from concurrent.futures import ProcessPoolExecutor
import pysam
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from celescope.tools import utils
from celescope.tools.step import Step, s_common
from celescope.trust_vdj.__init__ import INDEX, TOOLS_DIR, CHAIN, CONDA_PATH
from celescope.tools.cellranger3 import get_plot_elements


@utils.add_log
def trust_assemble(thread, index, fq, bc_fa, umi_fa, outdir, sample, trimLevel=1):

    cmd = (
        f'{CONDA_PATH}/bin/trust4  -t {thread} '
        f'-f {index} '
        f'-o {outdir}/{sample} '
        f'-u {fq} '
        f'--barcode {bc_fa} '
        f'--UMI {umi_fa} '
        f'--trimLevel {trimLevel}'
    )
    trust_assemble.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)

@utils.add_log
def get_full_len_assembly(outdir, sample, annotated_fa):
    cmd = (
        f'perl {TOOLS_DIR}/GetFullLengthAssembly.pl '
        f'{outdir}/{sample}_annot.fa > '
        f'{outdir}/{sample}_full_len.fa '
    )
    get_full_len_assembly.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def annotate(index, input_file, sample, thread, outdir, species):

    cmd = (
        f'{CONDA_PATH}/bin/annotator -f {INDEX}/{species}/IMGT+C.fa '
        f'-a {outdir}/{sample}_final.out '
        f'-t {thread} '
        f'-o {outdir}/{sample} '
        f'--barcode --UMI '
        f'--readAssignment {outdir}/{sample}_assign.out '
        f'-r {outdir}/{sample}_assembled_reads.fa > {outdir}/{sample}_annot.fa'
    )
    annotate.logger.info(cmd)
    subprocess.check_call(cmd, shell=True)


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
        self.temp_dir = f'{self.assemble_out}/temp'
        # check dir
        utils.check_mkdir(self.match_out)
        utils.check_mkdir(self.assemble_out)
        utils.check_mkdir(self.summarize_out)
        utils.check_mkdir(self.temp_dir)

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
    def process(self):
        # process all vdj
        cb_range = self.barcodeRange.split(' ')
        umi_range = self.umiRange.split(' ')
        cmd = (
            f'{CONDA_PATH}/bin/fastq-extractor -t {self.thread} '
            f'-f {INDEX}/{self.species}/bcrtcr.fa '
            f'-o {self.temp_dir}/{self.sample} '
            f'--barcodeStart {int(cb_range[0])} '
            f'--barcodeEnd {int(cb_range[1])} '
            f'--umiStart {int(umi_range[0])} '
            f'--umiEnd {int(umi_range[1])} '
            f'-u {self.match_fq2} '
            f'--barcode {self.match_fq1} '
            f'--UMI {self.match_fq1} ' 
        )
        Assemble.process.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

        # process single chain
        for c in self.chains:
            cmd = (
                f'{CONDA_PATH}/bin/fastq-extractor -t {self.thread} '
                f'-f {INDEX}/{self.species}/{c}.fa '
                f'-o {self.temp_dir}/{c}_toassemble '
                f'--barcodeStart {int(cb_range[0])} '
                f'--barcodeEnd {int(cb_range[1])} '
                f'--umiStart {int(umi_range[0])} '
                f'--umiEnd {int(umi_range[1])} '
                f'-u {self.match_fq2} '
                f'--barcode {self.match_fq1} '
                f'--UMI {self.match_fq1} '
            )
            Assemble.process.logger.info(cmd)
            subprocess.check_call(cmd, shell=True)

            f = pysam.FastxFile(f'{self.temp_dir}/{c}_toassemble.fq')
            # self.cell_summary.append({
            #     'item': f'Reads Mapped to {c}', 
            #     'count': len(list(f)), 
            #     'total_count': self.matched_reads
            # })
            del f

        # cutoff by umi num
        fl = pysam.FastxFile(f'{self.temp_dir}/{self.sample}.fq')
        # self.cell_summary.append({
        #     'item': f'Reads Mapped to Any V(D)J genes', 
        #     'count': len(list(fl)),
        #     'total_count': self.matched_reads
        # })
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
        df_UMI_cell = df[df.UMI >= UMI_min]
        df["mark"] = df["UMI"].apply(
            lambda x: "CB" if (x >= UMI_min) else "UB")
        
        df.to_csv(f'{self.assemble_out}/count.txt', sep='\t', index=False)
        self.add_data_item(chart=get_plot_elements.plot_barcode_rank(f'{self.assemble_out}/count.txt'))

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

        references = [f'{INDEX}/{self.species}/bcrtcr.fa'] * len(idx)
        imgt_refs = [f'{INDEX}/{self.species}/IMGT+C.fa'] * len(idx)
        fqs = [f'{self.temp_dir}/temp_{i}.fq' for i in range(len(idx)-1)]
        bc_fas = [f'{self.temp_dir}/temp_{i}_bc.fa' for i in range(len(idx)-1)]
        umi_fas = [f'{self.temp_dir}/temp_{i}_umi.fa' for i in range(len(idx)-1)]
        threads = [self.thread] * len(idx)
        temp_dirs = [self.temp_dir] * len(idx)
        species = [self.species] * len(idx)
        samples = [f'temp_{i}' for i in range(len(idx)-1)]
        ass_outs = [f'{self.temp_dir}/temp_{i}_final.out' for i in range(len(idx)-1)]
        annot_outs = [f'{self.temp_dir}/temp_{i}_annot.fa' for i in range(len(idx)-1)]
        ass = []
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(trust_assemble, threads, references, fqs, bc_fas, umi_fas, temp_dirs, samples):
                ass.append(res)
        annot = []
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(annotate, imgt_refs, ass_outs, samples, threads, temp_dirs, species):
                annot.append(res) 
        gfla = []
        with ProcessPoolExecutor(self.thread) as pool:
            for res in pool.map(get_full_len_assembly, temp_dirs, samples, annot_outs):
                ass.append(res)
        # for i in range(len(idx)-1):
        #     fq = f'{self.temp_dir}/temp_{i}.fq'
        #     bc_fa = f'{self.temp_dir}/temp_{i}_bc.fa'
        #     umi_fa = f'{self.temp_dir}/temp_{i}_umi.fa'
        #     trust_assemble(self.thread, reference, fq, bc_fa, umi_fa, 
        #                     self.temp_dir, f'temp_{i}')
        #     annotate(imgt_ref, f'{self.temp_dir}/temp_{i}_final.out', 
        #                 f'temp_{i}', self.thread, self.temp_dir, self.species)
        #     get_full_len_assembly(self.temp_dir, f'temp_{i}', f'{self.temp_dir}/temp_{i}_annot.fa')

        temp_outs_fa = glob.glob(f'{self.temp_dir}/temp_*_annot.fa')
        string = ' '.join(temp_outs_fa)
        cmd = f'cat {string} > {self.assemble_out}/{self.sample}_annot.fa'
        Assemble.process.logger.info(cmd)
        os.system(cmd)

    def run(self):
        if not os.path.exists(self.match_fq2):
            self.get_match_fastq()
        self.process()
        os.system(f'rm -rf {self.temp_dir}')

@utils.add_log
def assemble(args):
    step_name = 'assemble'
    assemble_obj = Assemble(args, step_name)
    assemble_obj.run()


def get_opts_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq2', help='R2 reads matched with scRNA-seq.', required=True)
        parser.add_argument('--match_dir', help='match scRNA-seq directory.', required=True)

    parser.add_argument('--species', help='species.', choices=["hg19", "hg38", "GRCm38"], required=True)
    parser.add_argument('--seqtype', help='TCR/BCR seq data.', choices=['TCR', 'BCR'], required=True)
    parser.add_argument('--barcodeRange', help='barcode range in fq1, INT INT CHAR.', default='0 23 +') 
    parser.add_argument('--umiRange', help='umi range in fq1, INT INT CHAR.', default='24 -1 +')
    parser.add_argument('--trimLevel', help='INT: 0: no trim; 1: trim low quality; 2: trim unmatched.', default=1)
    parser.add_argument('--expect_cells', help='Expected Cells number, INT.', default=3000)








