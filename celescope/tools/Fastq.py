from collections import defaultdict
import pysam
import gzip
import numpy as np
import subprocess
from xopen import xopen
from concurrent.futures import ProcessPoolExecutor
from celescope.tools.utils import genDict, log, fastq_line

class Fastq():
    """
    Manipulate fastq file
    """

    def __init__(self, fq):
        self.fq = fq
    
    @log
    def sort_fastq(self, fq_tmp_file):
        cmd = (
            f'zcat {self.fq} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {fq_tmp_file};'
        )
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    @log
    def sorted_dumb_consensus(fq, outfile):
        '''
        read in name sorted fastq, output (barcode,umi) consensus fastq
        '''
        curr_combine = ""
        read_list = []
        n = 0
        out_h = xopen(outfile, 'w')

        with pysam.FastxFile(fq) as fh:
            for entry in fh:
                attr = entry.name.split('_')
                barcode = attr[0]
                umi = attr[1]
                combine = [barcode,umi]
                if combine != curr_combine:
                    # first
                    if curr_combine == "":
                        curr_combine = combine
                        read_list.append([entry.sequence,entry.quality])
                        continue
                    consensus, consensus_qual = Fastq.dumb_consensus(read_list, threshold=0.5, ambiguous="N")
                    n += 1
                    prefix = "_".join(curr_combine)
                    read_name = f'{prefix}_{n}'
                    out_h.write(fastq_line(read_name,consensus,consensus_qual))
                    if n % 10000 == 0:
                        Fastq.sorted_dumb_consensus.logger.info(f'{n} UMI done.')
                    read_list = []
                    curr_combine = combine
                read_list.append([entry.sequence,entry.quality])
        #last
        consensus, consensus_qual = Fastq.dumb_consensus(read_list, threshold=0.5, ambiguous="N")  
        n += 1
        read_name = f'{barcode}_{umi}_{n}'
        out_h.write(fastq_line(read_name,consensus,consensus_qual))

        out_h.close()
        return n

    @log
    def wrap_consensus(self, outdir, sample):
        fq_tmp_file = f'{outdir}/{sample}_sorted.fq.tmp'
        self.sort_fastq(fq_tmp_file)
        outfile = f'{outdir}/{sample}_consensus.fq'
        n = self.sorted_dumb_consensus(fq=fq_tmp_file, outfile=outfile)
        return outfile, n


    @log
    def split_fq(self):
        '''
        split reads into barcode_UMI dict
        read_list = fq_dict[barcode][umi]
        elements of read_list: [entry.sequence,entry.quality]
        '''
        fq_dict = genDict(dim=2, valType=list)

        with pysam.FastxFile(self.fq) as fh:
            for entry in fh:
                attr = entry.name.split('_')
                barcode = attr[0]
                umi = attr[1]
                fq_dict[barcode][umi].append([entry.sequence,entry.quality])
        return fq_dict

    @staticmethod
    def dumb_consensus(read_list, threshold=0.5, ambiguous="N"):
        '''
        This is similar to biopython dumb_consensus.
        It will just go through the sequence residue by residue and count up the number of each type
        of residue (ie. A or G or T or C for DNA) in all sequences in the
        alignment. If the percentage of the most common residue type is
        greater then the passed threshold, then we will add that residue type,
        otherwise an ambiguous character will be added.
        '''

        con_len = Fastq.get_read_length(read_list, threshold=0.5)
        consensus = ""
        consensus_qual = ""
        for n in range(con_len):
            atom_dict = defaultdict(int)
            quality_dict = defaultdict(int)
            num_atoms = 0
            for read in read_list:
                # make sure we haven't run past the end of any sequences
                # if they are of different lengths
                sequence = read[0]
                quality = read[1]
                if n < len(sequence):
                    atom = sequence[n]
                    atom_dict[atom] += 1
                    num_atoms = num_atoms + 1

                    base_qual = quality[n]
                    quality_dict[base_qual] += 1

            consensus_atom = ambiguous
            for atom in atom_dict:
                if atom_dict[atom] >= num_atoms * threshold:
                    consensus_atom = atom
                    break
            consensus += consensus_atom

            max_freq_qual = 0
            consensus_base_qual = 'F'
            for base_qual in quality_dict:
                if quality_dict[base_qual] > max_freq_qual:
                    max_freq_qual = quality_dict[base_qual]
                    consensus_base_qual = base_qual

            consensus_qual += consensus_base_qual
        return consensus, consensus_qual

    @staticmethod
    def get_read_length(read_list, threshold=0.5):
        '''
        compute read_length from read_list. 
        length = max length with read fraction >= threshold
        '''
        
        n_read = len(read_list)
        length_dict = defaultdict(int)
        for read in read_list:
            length = len(read[0])
            length_dict[length] += 1
        for length in length_dict:
            length_dict[length] = length_dict[length] / n_read

        fraction = 0
        for length in sorted(length_dict.keys(),reverse=True):
            fraction += length_dict[length]
            if fraction >= threshold:
                return length
    
    @log
    def umi_dumb_consensus(self, threshold=0.5):
        '''
        dumb_consensus for barcode,umi
        '''
        consensus_dict = genDict(dim=2, valType = list)
        fq_dict = self.split_fq()
        n_umi = 0
        for barcode in fq_dict:
            for umi in fq_dict[barcode]:
                consenesus, consensus_qual = Fastq.dumb_consensus(fq_dict[barcode][umi], threshold=threshold)
                consensus_dict[barcode][umi] = [consenesus, consensus_qual]
                n_umi += 1
                if n_umi % 10000 == 0:
                    Fastq.umi_dumb_consensus.logger.info(f'{n_umi} UMI done.')
        self.consensus_dict = consensus_dict

    @staticmethod
    @log
    def dumb_consensus_worker(fq_dict, threshold):
        consensus_dict = {}        
        for barcode in fq_dict:
            for umi in fq_dict[barcode]:
                consenesus, consensus_qual = Fastq.dumb_consensus(fq_dict[barcode][umi], threshold=threshold)
                if barcode not in consensus_dict:
                    consensus_dict[barcode] = {}
                consensus_dict[barcode][umi] = [consenesus, consensus_qual]
        return consensus_dict

    @staticmethod
    @log 
    def split_dict(input_dict: dict, num_parts: int) -> list:
        return [dict(list(input_dict.items())[i::num_parts])
            for i in range(num_parts)]


    @log
    def umi_dumb_consensus_concurrent(self, thread, threshold=0.5):
        '''
        dumb_consensus for barcode,umi
        '''

        fq_dict = self.split_fq()
        consensus_dict = {}

        fq_dicts = Fastq.split_dict(fq_dict, thread)
        for fq_dict in fq_dicts:
            umi_list = []
            read_list = []
            for barcode in fq_dict.keys():
                umi_list.append(len(fq_dict[barcode]))
                for umi in fq_dict[barcode]:
                    read_list.append(len(fq_dict[barcode][umi]))

            print('mean UMI count:' + str(np.mean(umi_list)))
            print('mean read count:' +  str(np.mean(read_list)))

        threshold_list = [threshold] * thread

        with ProcessPoolExecutor(thread) as pool:
            for res in pool.map(Fastq.dumb_consensus_worker, fq_dicts, threshold_list):
                consensus_dict.update(res)
        self.consensus_dict = consensus_dict

    @log
    def write_consensus_fasta(self, outdir, sample):
        out_fasta = f'{outdir}/{sample}_consensus.fasta'
        index = 0
        import os
        print(os.getcwd())
        with open(out_fasta, 'wt') as handle:
            for barcode in self.consensus_dict:
                for umi in self.consensus_dict[barcode]:
                    index += 1
                    read_name = f'{barcode}_{umi}_{index}'
                    seq = self.consensus_dict[barcode][umi][0]
                    handle.write('>' + read_name + '\n')
                    handle.write(seq + '\n')
        return out_fasta

    @log
    def write_consensus_fastq(self, outdir, sample, use_gzip=True):
        suffix = ''
        if gzip:
            suffix = '.gz'
        out_fastq = f'{outdir}/{sample}_consensus.fastq{suffix}'
        index = 0
        with gzip.open(out_fastq, 'wt') as handle:
            for barcode in self.consensus_dict:
                for umi in self.consensus_dict[barcode]:
                    index += 1
                    read_name = f'{barcode}_{umi}_{index}'
                    seq = self.consensus_dict[barcode][umi][0]
                    qual = self.consensus_dict[barcode][umi][1]
                    handle.write('@' + read_name + '\n')
                    handle.write(seq + '\n')
                    handle.write('+\n')
                    handle.write(qual + '\n')
        return out_fastq



                

        