#!/bin/env python
# coding=utf8

import os
import re
import sys
import json
import logging
import subprocess
import glob
from celescope.tools.utils import format_number
from celescope.tools.utils import glob_genomeDir
from celescope.tools.report import reporter

logger1 = logging.getLogger(__name__)


def format_stat(map_log, region_log, samplename):
    fh1 = open(map_log, 'r')
    p1 = re.compile(r'Uniquely mapped reads number\s+(\d+)')
    UNIQUE_READS = []
    MULTI_MAPPING_READS = []
    for line in fh1:
        if line.strip() == '':
            continue
        if re.search(r'Uniquely mapped reads', line):
            UNIQUE_READS.append(line.strip().split()[-1])
        if re.search(r'of reads mapped to too many loci', line):
            MULTI_MAPPING_READS.append(line.strip().split()[-1])

    fh2 = open(region_log, 'r')
    region_dict = {}
    while True:
        line = fh2.readline()
        if not line:
            break
        if line.startswith('## METRICS CLASS'):
            header = fh2.readline().strip().split('\t')
            data = fh2.readline().strip().split('\t')
            region_dict = dict(zip(header, data))
            break

    Total = float(region_dict['PF_ALIGNED_BASES'])
    Exonic_Regions = int(region_dict['UTR_BASES']) + \
        int(region_dict['CODING_BASES'])
    Intronic_Regions = int(region_dict['INTRONIC_BASES'])
    Intergenic_Regions = int(region_dict['INTERGENIC_BASES'])

    region_dict['Exonic_Regions'] = "{}({:.2%})".format(
        format_number(Exonic_Regions), Exonic_Regions / Total)
    region_dict['Intronic_Regions'] = "{}({:.2%})".format(
        format_number(Intronic_Regions), Intronic_Regions / Total)
    region_dict['Intergenic_Regions'] = "{}({:.2%})".format(
        format_number(Intergenic_Regions), Intergenic_Regions / Total)

    with open(os.path.dirname(map_log) + '/stat.txt', 'w') as stat_fh:
        stat_fh.write('%s: %s(%s)\n' % ('Uniquely Mapped Reads',
                                        format_number(int(UNIQUE_READS[0])), UNIQUE_READS[1]))
        stat_fh.write('%s: %s(%s)\n' % ('Multi-Mapped Reads',
                                        format_number(
                                            int(MULTI_MAPPING_READS[0])),
                                        MULTI_MAPPING_READS[1]))
        stat_fh.write(
            '%s: %s\n' %
            ('Base Pairs Mapped to Exonic Regions',
             region_dict['Exonic_Regions']))
        stat_fh.write(
            '%s: %s\n' %
            ('Base Pairs Mapped to Intronic Regions',
             region_dict['Intronic_Regions']))
        stat_fh.write(
            '%s: %s\n' %
            ('Base Pairs Mapped to Intergenic Regions',
             region_dict['Intergenic_Regions']))
    return {'region_labels': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
            'region_values': [Exonic_Regions, Intronic_Regions, Intergenic_Regions]}


def STAR(args):
    logger1.info('STAR ...!')
    # check
    refFlat, gtf = glob_genomeDir(args.genomeDir, logger1)

    # check dir
    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % (args.outdir))

    # run STAR
    outPrefix = args.outdir + '/' + args.sample + '_'
    outBam = args.outdir + '/' + args.sample + '_'
    # cmd = ['STAR', '--runThreadN', str(args.thread), '--genomeDir', args.genomeDir, '--readFilesIn', args.fq, '--readFilesCommand', 'zcat', '--outFilterMultimapNmax', '1', '--outReadsUnmapped', 'Fastx', '--outFileNamePrefix', outPrefix, '--outSAMtype', 'BAM', 'SortedByCoordinate']
    cmd = ['STAR', '--runThreadN', str(args.thread), '--genomeDir', args.genomeDir,
           '--readFilesIn', args.fq, '--readFilesCommand', 'zcat', '--outFilterMultimapNmax',
           '1', '--outFileNamePrefix', outPrefix, '--outSAMtype', 'BAM', 'SortedByCoordinate']
    if args.out_unmapped:
        cmd.append(['--outReadsUnmapped', 'Fastx'])
    logger1.info('%s' % (' '.join(cmd)))
    subprocess.check_call(cmd)
    logger1.info('STAR done!')

    logger1.info('stat mapping region ...!')
    outBam = outPrefix + 'Aligned.sortedByCoord.out.bam'
    region_txt = args.outdir + '/' + args.sample + '_region.log'
    cmd = [
        'picard',
        '-Xmx4G',
        '-XX:ParallelGCThreads=4',
        'CollectRnaSeqMetrics',
        'I=%s' %
        (outBam),
        'O=%s' %
        (region_txt),
        'REF_FLAT=%s' %
        (refFlat),
        'STRAND=NONE',
        'VALIDATION_STRINGENCY=SILENT']
    logger1.info('%s' % (' '.join(cmd)))
    res = subprocess.run(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    logger1.info(res.stdout)
    logger1.info('stat mapping region done!')

    plot = format_stat(
        args.outdir +
        '/' +
        args.sample +
        '_Log.final.out',
        region_txt,
        args.sample)
    t = reporter(
        name='STAR',
        assay=args.assay,
        sample=args.sample,
        stat_file=args.outdir + '/stat.txt',
        outdir=args.outdir + '/..',
        plot=plot)
    t.get_report()


def get_opts_STAR(parser, sub_program):
    if sub_program:
        parser.add_argument('--fq', required=True)
        parser.add_argument('--readFilesCommand', default='zcat')
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument('--thread', default=1)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--out_unmapped',
        help='out_unmapped',
        action='store_true')
    parser.add_argument('--genomeDir', help='genome directory', required=True)
