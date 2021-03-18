from collections import defaultdict
import pysam
import gzip
import numpy as np
import subprocess
import os
from xopen import xopen
from celescope.tools.utils import format_metrics, format_ratios, log, fastq_line, gen_stat
from celescope.tools.report import reporter


@log
def sort_fastq(fq, fq_tmp_file, outdir):
    tmp_dir = f'{outdir}/tmp'
    cmd = (
        f'mkdir {tmp_dir};'
        f'zcat {fq} | paste - - - - | sort -T {tmp_dir} -k1,1 -t " " | tr "\t" "\n" > {fq_tmp_file};'
    )
    subprocess.check_call(cmd, shell=True)


@log
def sorted_dumb_consensus(fq, outfile, threshold):
    '''
    read in name sorted fastq, output (barcode,umi) consensus fastq
    '''
    curr_combine = ""
    read_list = []
    n = 0
    total_ambiguous_base_n = 0
    length_list = []
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
                consensus, consensus_qual, ambiguous_base_n, con_len = dumb_consensus(read_list, threshold=threshold, ambiguous="N")
                n += 1
                prefix = "_".join(curr_combine)
                read_name = f'{prefix}_{n}'
                out_h.write(fastq_line(read_name,consensus,consensus_qual))
                if n % 10000 == 0:
                    sorted_dumb_consensus.logger.info(f'{n} UMI done.')
                total_ambiguous_base_n += ambiguous_base_n
                length_list.append(con_len)
                read_list = []
                curr_combine = combine
            read_list.append([entry.sequence,entry.quality])
    #last
    consensus, consensus_qual, ambiguous_base_n, con_len = dumb_consensus(read_list, threshold=threshold, ambiguous="N")
    n += 1
    prefix = "_".join(curr_combine)
    read_name = f'{prefix}_{n}'
    out_h.write(fastq_line(read_name,consensus,consensus_qual))
    if n % 10000 == 0:
        sorted_dumb_consensus.logger.info(f'{n} UMI done.')
    total_ambiguous_base_n += ambiguous_base_n
    length_list.append(con_len)
    
    out_h.close()
    return n, total_ambiguous_base_n, length_list


@log
def wrap_consensus(fq, outdir, sample, threshold):
    fq_tmp_file = f'{outdir}/{sample}_sorted.fq.tmp'
    sort_fastq(fq, fq_tmp_file, outdir)
    outfile = f'{outdir}/{sample}_consensus.fq'
    n, total_ambiguous_base_n, length_list = sorted_dumb_consensus(fq=fq_tmp_file, outfile=outfile, threshold=threshold)
    return outfile, n, total_ambiguous_base_n, length_list


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
    consensus = ""
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
        consensus += consensus_atom

        max_freq_qual = 0
        consensus_base_qual = default_qual
        for base_qual in quality_dict:
            if quality_dict[base_qual] > max_freq_qual:
                max_freq_qual = quality_dict[base_qual]
                consensus_base_qual = base_qual

        consensus_qual += consensus_base_qual
    return consensus, consensus_qual, ambiguous_base_n, con_len


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
    for length in sorted(length_dict.keys(),reverse=True):
        fraction += length_dict[length]
        if fraction >= threshold:
            return length


def consensus(args):
    sample = args.sample
    outdir = args.outdir
    assay = args.assay
    fq = args.fq
    threshold = float(args.threshold)

    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    outfile, n, total_ambiguous_base_n, length_list = wrap_consensus(fq, outdir, sample, threshold)

    # metrics
    metrics = {}
    metrics["UMI Counts"] = n
    metrics["Median UMI Length"] = np.median(length_list)
    metrics["Ambiguous Base Counts"] = total_ambiguous_base_n    
    format_metrics(metrics)

    ratios = {}
    ratios["Ambiguous Base Counts Ratio"] = total_ambiguous_base_n / sum(length_list)
    format_ratios(ratios)

    # stat file
    stat_file = f'{outdir}/stat.txt'
    with open(stat_file, 'w') as stat_h:
        stat_str = (
            f'UMI Counts: {metrics["UMI Counts"]}\n'
            f'Median UMI Length: {metrics["Median UMI Length"]}\n'
            f'Ambiguous Base Counts: {metrics["Ambiguous Base Counts"]}({ratios["Ambiguous Base Counts Ratio"]}%)\n'
        )
        stat_h.write(stat_str)

    t = reporter(name='consensus', assay=args.assay, sample=args.sample,
                 stat_file=stat_file, outdir=args.outdir + '/..')
    t.get_report()


def get_opts_consensus(parser, sub_program):
    parser.add_argument("--threshold", help='valid base threshold', default=0.5)
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--fq", required=True)
        parser.add_argument('--assay', help='assay', required=True)