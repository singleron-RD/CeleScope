import subprocess
from collections import defaultdict
from itertools import groupby

import numpy as np
import pysam
from xopen import xopen

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common


class Consensus(Step):
    """
    Features
    - Consensus all the reads of the same (barcode, UMI) combinations into one read(UMI).

    Output
    - `{sample}_consensus.fq` Consensus fastq.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # out files
        self.fq_tmp_file = f'{self.out_prefix}_sorted.fq.tmp'
        self.consensus_fq = f'{self.out_prefix}_consensus.fq'

    @utils.add_log
    def run(self):
        if self.args.not_consensus:
            Consensus.run.logger.warning("Will not perform UMI consensus!")
            return

        sort_fastq(self.args.fq, self.fq_tmp_file, self.outdir)
        n, total_ambiguous_base_n, length_list = sorted_dumb_consensus(
            fq=self.fq_tmp_file,
            outfile=self.consensus_fq,
            threshold=self.args.threshold
        )

        self.add_metric(
            name="UMI Counts",
            value=n,
        )
        self.add_metric(
            name="Mean UMI Length",
            value=round(np.mean(length_list), 2),
        )
        self.add_metric(
            name="Ambiguous Base Counts",
            value=total_ambiguous_base_n,
            total=sum(length_list),
        )
        self.clean_up()


@utils.add_log
def sort_fastq(fq, fq_tmp_file, outdir):
    tmp_dir = f'{outdir}/tmp'
    cmd = (
        f'mkdir {tmp_dir};'
        f'less {fq} | paste - - - - | sort -T {tmp_dir} -k1,1 -t " " | tr "\t" "\n" > {fq_tmp_file};'
    )
    subprocess.check_call(cmd, shell=True)


@utils.add_log
def sorted_dumb_consensus(fq, outfile, threshold):
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
                read_list, threshold=threshold, ambiguous="N")
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


def dumb_consensus(read_list, threshold=0.5, ambiguous='N', default_qual='F'):
    '''
    This is similar to biopython dumb_consensus.
    It will just go through the sequence residue by residue and count up the number of each type
    of residue (ie. A or G or T or C for DNA) in all sequences in the
    alignment. If the percentage of the most common residue type is
    greater then the passed threshold, then we will add that residue type,
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
            if atom_dict[atom] >= num_atoms * threshold:
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

    step_name = "consensus"
    consensus_obj = Consensus(args, step_name)
    consensus_obj.run()


def get_opts_consensus(parser, sub_program):
    parser.add_argument("--threshold", help='Default 0.5. Valid base threshold. ', type=float, default=0.5)
    parser.add_argument("--not_consensus", help="Skip the consensus step. ", action='store_true')
    if sub_program:
        parser.add_argument("--fq", help="Required. Fastq file.", required=True)
        s_common(parser)
