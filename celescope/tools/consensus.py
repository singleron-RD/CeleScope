import subprocess
import unittest
from collections import defaultdict
from itertools import groupby

import numpy as np
import pysam
from xopen import xopen

from celescope.tools import utils
from celescope.tools.step import Step, s_common


class Consensus(Step):
    """
    ## Features
    - Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI). It will go through the sequence residue by residue and 
    count up the number of each type of residue (ie. A or G or T or C for DNA) in all sequences in the
    alignment. If the following conditions are met, the consensus sequence will be the most common residue in the alignment:
    1. the percentage of the most common residue type > threshold(default: 0.5);
    2. most common residue reads >= min_consensus_read;
    otherwise an ambiguous character(N) will be added.

    ## Output
    - `{sample}_consensus.fq` Fastq file after consensus.
    """

    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        # set
        self.min_consensus_read = int(self.args.min_consensus_read)

        # out files
        self.fq_tmp_file = f'{self.out_prefix}_sorted.fq.tmp'
        self.consensus_fq = f'{self.out_prefix}_consensus.fq'

    @utils.add_log
    def run(self):

        sort_fastq(self.args.fq, self.fq_tmp_file, self.outdir)
        n, total_ambiguous_base_n, length_list = sorted_dumb_consensus(
            fq=self.fq_tmp_file,
            outfile=self.consensus_fq,
            threshold=self.args.threshold,
            min_consensus_read=self.min_consensus_read,
        )

        self.add_metric(
            name="UMI Counts",
            value=n,
            help_info='total UMI from FASTQ files',
        )
        self.add_metric(
            name="Mean UMI Length",
            value=round(np.mean(length_list), 2),
            help_info='mean of all UMI length'
        )
        self.add_metric(
            name="Ambiguous Base Counts",
            value=total_ambiguous_base_n,
            total=sum(length_list),
            help_info='number of bases that do not pass consensus threshold'
        )


@utils.add_log
def sort_fastq(fq, fq_tmp_file, outdir):
    tmp_dir = f'{outdir}/tmp'
    cmd = (
        f'mkdir {tmp_dir};'
        f'less {fq} | paste - - - - | sort -T {tmp_dir} -k1,1 -t " " | tr "\t" "\n" > {fq_tmp_file};'
    )
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def sorted_dumb_consensus(fq, outfile, threshold, min_consensus_read):
    '''
    consensus read in name-sorted fastq
    output (barcode,umi) consensus fastq
    '''
    read_list = []
    n_umi = 0
    total_ambiguous_base_n = 0
    length_list = []
    out_h = xopen(outfile, 'w')

    def keyfunc(read):
        attr = read.name.split('_')
        return (attr[0], attr[1])

    with pysam.FastxFile(fq) as fh:
        for (barcode, umi), g in groupby(fh, key=keyfunc):
            read_list = []
            for read in g:
                read_list.append([read.sequence, read.quality])
            consensus_seq, consensus_qual, ambiguous_base_n, con_len = dumb_consensus(
                read_list,
                threshold=threshold,
                min_consensus_read=min_consensus_read,
                ambiguous="N"
            )
            n_umi += 1
            prefix = "_".join([barcode, umi])
            read_name = f'{prefix}_{n_umi}'
            out_h.write(utils.fastq_line(read_name, consensus_seq, consensus_qual))
            if n_umi % 10000 == 0:
                sorted_dumb_consensus.logger.info(f'{n_umi} UMI done.')
            total_ambiguous_base_n += ambiguous_base_n
            length_list.append(con_len)

    out_h.close()
    return n_umi, total_ambiguous_base_n, length_list


def dumb_consensus(read_list, threshold=0.5, min_consensus_read=1, ambiguous='N', default_qual='F'):
    '''
    This is similar to biopython dumb_consensus.
    It will just go through the sequence residue by residue and count up the number of each type
    of residue (ie. A or G or T or C for DNA) in all sequences in the
    alignment. If 
    1. the percentage of the most common residue type > threshold;
    2. most common residue reads >= min_consensus_read;
    then we will add that residue type,
    otherwise an ambiguous character will be added.
    elements of read_list: [entry.sequence,entry.quality]
    '''

    con_len = get_read_length(read_list, threshold=threshold)
    consensus_seq = ""
    consensus_qual = ""
    ambiguous_base_n = 0
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
            if atom_dict[atom] > num_atoms * threshold and atom_dict[atom] >= min_consensus_read:
                consensus_atom = atom
                break
        if consensus_atom == ambiguous:
            ambiguous_base_n += 1
        consensus_seq += consensus_atom

        max_freq_qual = 0
        consensus_base_qual = default_qual
        for base_qual in quality_dict:
            if quality_dict[base_qual] > max_freq_qual:
                max_freq_qual = quality_dict[base_qual]
                consensus_base_qual = base_qual

        consensus_qual += consensus_base_qual
    return consensus_seq, consensus_qual, ambiguous_base_n, con_len


def get_read_length(read_list, threshold=0.5):
    '''
    compute read_length from read_list. 
    length = max length with read fraction >= threshold
    elements of read_list: [entry.sequence,entry.quality]
    '''

    n_read = len(read_list)
    length_dict = defaultdict(int)
    for read in read_list:
        length = len(read[0])
        length_dict[length] += 1
    for length in length_dict:
        length_dict[length] = length_dict[length] / n_read

    fraction = 0
    for length in sorted(length_dict.keys(), reverse=True):
        fraction += length_dict[length]
        if fraction >= threshold:
            return length


@utils.add_log
def consensus(args):
    if args.not_consensus:
        return
    with Consensus(args, display_title="Consensus") as runner:
        runner.run()


def get_opts_consensus(parser, sub_program):
    parser.add_argument("--threshold", help='Default 0.5. Valid base threshold. ', type=float, default=0.5)
    parser.add_argument("--not_consensus", help="Skip the consensus step. ", action='store_true')
    parser.add_argument("--min_consensus_read", help="Minimum number of reads to support a base. ", default=1)
    if sub_program:
        parser.add_argument("--fq", help="Required. Fastq file.", required=True)
        s_common(parser)


class Consensus_unittest(unittest.TestCase):

    def test_get_read_length(self):
        read_list = [['AAAA', 'FFFF'], ['TTT', 'FFF'], ['CCC', 'FFF'], ['GGGGGGG', 'FFFFFFF']]
        assert get_read_length(read_list, 0.5) == 4

    def test_dumb_consensus(self):
        read_list = [('AAAA', 'FFFF'), ('TTT', 'FF;'), ('CCCC', 'FFFF'), ('GGGAGGG', 'FFFFFFF')]
        consensus_seq, _consensus_qual, _ambiguous_base_n, _con_len = dumb_consensus(read_list, 0.5)
        self.assertEqual(consensus_seq, 'NNNA')

        read_list = [('AAAA', 'FFFF'), ('TTT', 'FF;'), ('CCC', 'FFF'), ('GGGGGGG', 'FFFFFFF')]
        consensus_seq, _consensus_qual, _ambiguous_base_n, _con_len = dumb_consensus(read_list, 0.5)
        self.assertEqual(consensus_seq, 'NNNN')


if __name__ == '__main__':
    unittest.main()
