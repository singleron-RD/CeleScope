#!/bin/env python
# coding=utf8

import os
import logging

logger1 = logging.getLogger(__name__)


def STAR_fusion(args):

    logger1.info('STAR_fusion...')

    sample = args.sample
    outdir = args.outdir
    input_read = args.input_read
    genomeDir = args.genomeDir
    runThreadN = args.thread

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % (outdir))

    out_prefix = f'{outdir}/{sample}_'
    out_BAM = out_prefix + "Aligned.sortedByCoord.out.bam"

    cmd = f"STAR \
 --genomeDir {genomeDir} \
 --readFilesIn {input_read}\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --runThreadN {runThreadN}\
 --limitBAMsortRAM 10000000000\
 --outFileNamePrefix {out_prefix}"

    logger1.info(cmd)
    os.system(cmd)
    logger1.info("STAR fusion done.")

    cmd = "samtools index {out_BAM}".format(out_BAM=out_BAM)
    logger1.info(cmd)
    os.system(cmd)
    logger1.info("samtools index done.")


def get_opts_STAR_fusion(parser, sub_program):
    if sub_program:
        parser.add_argument('--outdir', help='output dir', required=True)
        parser.add_argument('--sample', help='sample name', required=True)
        parser.add_argument("--input_read", required=True)
        parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument(
        '--genomeDir',
        help='fusion gene STAR index genome directory',
        required=True)
    parser.add_argument("--thread", help='STAR thread', default=1)
