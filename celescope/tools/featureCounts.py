#!/bin/env python
# coding=utf8

import os
import re
import subprocess
import pysam
from celescope.tools.utils import add_log, format_number, gene_convert, glob_genomeDir, s_common
from celescope.tools.step import Step


def format_stat(log):
    #Assigned, Unassigned_NoFeatures, Unassigned_Ambiguity=(0, 0, 0)
    tmp_arr = []
    fh = open(log, 'r')
    with open(os.path.dirname(log) + '/stat.txt', 'w') as stat_fh:
        p1 = re.compile(r'Assigned.*?(\d+)', flags=re.S)
        p2 = re.compile(r'Unassigned_NoFeatures.*?(\d+)', flags=re.S)
        p3 = re.compile(r'Unassigned_Ambiguity.*?(\d+)', flags=re.S)
        for line in fh:
            if line.strip() == '':
                continue

            m1 = p1.search(line.strip())
            if m1:
                tmp_arr.append(int(m1.group(1)))

            m2 = p2.search(line)
            if m2:
                tmp_arr.append(int(m2.group(1)))

            m3 = p3.search(line)
            if m3:
                tmp_arr.append(int(m3.group(1)))

        total = sum(tmp_arr)
        tmp_arr = [
            '%s(%.2f%%)' %
            (format_number(n), (n + 0.0) / total * 100) for n in tmp_arr]
        for t, s in zip(['Assigned', 'Unassigned_NoFeatures',
                         'Unassigned_Ambiguity'], tmp_arr):
            stat_fh.write('%s: %s\n' % (t, s))
    fh.close()

@add_log
def add_tag(bam, gtf):
    id_name = gene_convert(gtf)
    samfile = pysam.AlignmentFile(bam, "rb")
    header = samfile.header
    new_bam = pysam.AlignmentFile(
        bam + ".temp", "wb", header=header)
    for read in samfile:
        attr = read.query_name.split('_')
        barcode = attr[0]
        umi = attr[1]
        read.set_tag(tag='CB', value=barcode, value_type='Z')
        read.set_tag(tag='UB', value=umi, value_type='Z')
        if read.has_tag('XT'):
            gene_id = read.get_tag('XT')
            gene_name = id_name[gene_id]
            read.set_tag(tag='GN', value=gene_name, value_type='Z')
            read.set_tag(tag='GX', value=gene_id, value_type='Z')
        new_bam.write(read)
    new_bam.close()
    cmd = f'mv {bam}.temp {bam}'
    subprocess.check_call(cmd, shell=True)

@add_log
def featureCounts(args):

    step_name = "featureCounts"
    step = Step(args, step_name)

    # check
    if args.genomeDir and args.genomeDir != "None":
        _refFlat, gtf, _ = glob_genomeDir(args.genomeDir)
    else:
        gtf = args.gtf

    # run featureCounts
    outPrefix = args.outdir + '/' + args.sample
    cmd = [
        'featureCounts',
        '-s', '1',
        '-a', gtf,
        '-o', outPrefix,
        '-R', 'BAM',
        '-T', str(args.thread),
        '-t', args.gtf_type,
        args.input,
    ]
    featureCounts.logger.info('%s' % (' '.join(cmd)))
    subprocess.check_call(cmd)

    subprocess.check_call(['which', 'samtools'])

    bam_basename = os.path.basename(args.input)
    bam = f'{args.outdir}/{bam_basename}.featureCounts.bam'

    # add tag
    add_tag(bam, gtf)

    # sort by name:BC and umi
    featureCounts.logger.info('samtools sort by name...!')
    cmd = [
        'samtools',
        'sort',
        '-n',
        '-@',
        str(args.thread),
        '-o',
        outPrefix +
        '_name_sorted.bam',
        bam]
    featureCounts.logger.info('%s' % (' '.join(cmd)))
    subprocess.check_call(cmd)
    featureCounts.logger.info('samtools sort sort by name done.')

    format_stat(args.outdir + '/' + args.sample + '.summary')

    step.clean_up()


def get_opts_featureCounts(parser, sub_program):
    parser.add_argument('--gtf', help='gtf file path')
    parser.add_argument(
        '--gtf_type',
        help='Specify feature type in GTF annotation',
        default='exon')
    parser.add_argument('--genomeDir')
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--input', required=True)
    return parser

