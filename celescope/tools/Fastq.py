from collections import defaultdict
import pysam
import gzip
from celescope.tools.utils import genDict, log

class Fastq():
    """
    Manipulate fastq from 02.cutadapt
    """

    def __init__(self, fq):
        self.fq = fq

    @log
    def split_fq(self):
        '''
        split barcode into barcode_UMI dict
        '''
        fq_dict = genDict(dim=2, valType=list)

        with pysam.FastxFile(self.fq) as fh:
            for entry in fh:
                attr = entry.name.split('_')
                barcode = attr[0]
                umi = attr[1]
                fq_dict[barcode][umi].append(entry.sequence)
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
        for n in range(con_len):
            atom_dict = defaultdict(int)
            num_atoms = 0
            for read in read_list:
                # make sure we haven't run past the end of any sequences
                # if they are of different lengths
                if n < len(read):
                    atom = read[n]
                    atom_dict[atom] += 1
                    num_atoms = num_atoms + 1

            consensus_atom = ambiguous
            for atom in atom_dict:
                if atom_dict[atom] >= num_atoms * threshold:
                    consensus_atom = atom
                    break 
            consensus += consensus_atom
        return consensus

    @staticmethod
    def get_read_length(read_list, threshold=0.5):
        '''
        compute read_length from read_list. 
        length = max length with read fraction >= threshold 

        '''
        n_read = len(read_list)
        length_dict = defaultdict(int)
        for read in read_list:
            length = len(read)
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
        consensus_dict = genDict(dim=2, valType = str)
        fq_dict = self.split_fq()
        for barcode in fq_dict:
            for umi in fq_dict[barcode]:
                consenesus = Fastq.dumb_consensus(fq_dict[barcode][umi], threshold=threshold)
                consensus_dict[barcode][umi] = consenesus
        self.consensus_dict = consensus_dict

    @log
    def write_consensus_fasta(self, outdir, sample):
        out_fasta = f'{outdir}/{sample}_consensus.fasta.gz'
        index = 0
        with gzip.open(out_fasta, 'wt') as handle:
            for barcode in self.consensus_dict:
                for umi in self.consensus_dict[barcode]:
                    index += 1
                    read_name = f'{barcode}_{umi}_{index}'
                    seq = self.consensus_dict[barcode][umi]
                    handle.write(read_name + '\n')
                    handle.write(seq + '\n')



                

        