#!/bin/env python
#coding=utf8

import os, re
import logging
import subprocess

FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(level = logging.INFO, format = FORMAT)

def get_opts4(parser,sub_program): 

    parser.add_argument('--thread', default=2)
    parser.add_argument('--annot', required=True)
    if sub_program:
        parser.add_argument('--input', required=True)
        #parser.add_argument('--format', default='BAM')
        parser.add_argument('--outdir', help='output dir',required=True)
        parser.add_argument('--sample', help='sample name', required=True)


def format_stat(log, samplename):
    #Assigned, Unassigned_NoFeatures, Unassigned_Ambiguity=(0, 0, 0)
    tmp_arr = []
    fh = open(log, 'r')
    with open(os.path.dirname(log) + '/stat.txt', 'w') as stat_fh:
        stat_fh.write('%s: %s\n'%('SampleName', samplename))
        p1 = re.compile(r'Assigned.*?(\d+)', flags=re.S)
        p2 = re.compile(r'Unassigned_NoFeatures.*?(\d+)', flags=re.S)
        p3 = re.compile(r'Unassigned_Ambiguity.*?(\d+)', flags=re.S)
        for line in fh:
            if line.strip()=='': continue

            m1=p1.search(line.strip())
            if m1: tmp_arr.append(int(m1.group(1)))

            m2=p2.search(line)
            if m2: tmp_arr.append(int(m2.group(1)))

            m3=p3.search(line)
            if m3: tmp_arr.append(int(m3.group(1)))

        total = sum(tmp_arr)
        tmp_arr = ['%s(%.2f%%)'%(n, (n+0.0)/total*100) for n in tmp_arr]
        for t, s in zip(['Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity'], tmp_arr):
            stat_fh.write('%s: %s\n'%(t, s))
    fh.close()

def featureCounts(args):
    """
    """
    logging.info('featureCounts ...!')
    # check dir
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    # run featureCounts
    outPrefix = args.outdir + '/' + args.sample
    cmd = ['featureCounts','-t','gene','-a', args.annot, '-o', outPrefix, '-R', 'BAM', '-T', str(args.thread), args.input]
    logging.info('%s'%(' '.join(cmd)))
    subprocess.check_call(cmd)
    logging.info('featureCounts done!')

    subprocess.check_call(['which', 'samtools'])

    # sort by name:BC and umi 
    logging.info('samtools sort ...!')
    bam_basename = os.path.basename(args.input)
    cmd = ['samtools', 'sort', '-n', '-@','3', '-o', outPrefix+'_name_sorted.bam', args.outdir + '/' + bam_basename + '.featureCounts.bam']
    logging.info('%s'%(' '.join(cmd)))
    subprocess.check_call(cmd)
    logging.info('samtools sort done!')

    logging.info('generate report ...!')
    format_stat(args.outdir+'/'+args.sample+'.summary', args.sample)
    from report import reporter
    t = reporter(name='featureCounts', stat_file=args.outdir + '/stat.txt', outdir=args.outdir + '/..')
    t.get_report()
    logging.info('generate report done!')

def main():
    import argparse
    parser = argparse.ArgumentParser(description='')
    get_opts4(parser)
    args = parser.parse_args() 
    featureCounts(args)


if __name__ == '__main__':
    main()

