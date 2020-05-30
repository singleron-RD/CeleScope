#!/bin/env python
#coding=utf8

import os
import glob
import sys
import argparse
import re
from collections import defaultdict

toolsdir = os.path.realpath(sys.path[0] + '/../tools')

'''
def parse_map(mapfile):
    dict = defaultdict(list)
    with open(mapfile) as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            tmp = line.split('\t')
            dict[tmp[0]] = tmp[1:]

    return dict
'''

def parse_map(mapfile, cells=3000):
    fq_dict = defaultdict(list)
    cells_dict = defaultdict(list)
    sample_arr = []
    with open(mapfile) as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            tmp = line.split()
            try:
                fq1 = glob.glob(tmp[1] + '/' + tmp[0] + '*' + '_1.fq.gz')[0]
                fq2 = glob.glob(tmp[1] + '/' + tmp[0] + '*' + '_2.fq.gz')[0]
            except IndexError:
                sys.exit('glob err with %s'%(tmp[1] + '/' + tmp[0] + '_*' + '_[12].fq.gz'))
                
            assert os.path.exists(fq1), '%s not exists!'%(fq1)
            assert os.path.exists(fq2), '%s not exists!'%(fq2)
            fq_dict[tmp[2]] = [fq1, fq2]

            if re.match(r'\d+$', tmp[-1]):
                cells_dict[tmp[2]] = tmp[-1]
            else:
                cells_dict[tmp[2]] = cells
            sample_arr.append(tmp[2])

    return fq_dict, sample_arr, cells_dict

def generate_sjm(cmd, name, q='all.q', m=1, x=1):
    cmd = '''
job_begin
    name {name}
    sched_options -w n -cwd -V -l vf={m}g,p={x} -q {q}
    cmd {cmd}
job_end
'''.format(
    name = name, m=m, x=x, q=q, cmd=re.sub(r'\s+', r' ', cmd.replace('\n',' ')))

    return cmd

def main():
    parser = argparse.ArgumentParser('scope-tools for multisample')
    #parser.add_argument('--mod', help='mod, sjm or shell', choices=['sjm', 'shell'], default='sjm')
    parser.add_argument('--mapfile', help='mapfile, 3 columns, "LibName\\tDataDir\\tSampleName"', required=True)
    parser.add_argument('--whitelist', help='cellbarcode list')
    parser.add_argument('--linker', help='linker')
    parser.add_argument('--pattern', help='read1 pattern, default=C8L16C8L16C8U8T18', default='C8L16C8L16C8U8T18')
    parser.add_argument('--bcType', help='choice of barcode types. Currently support scope and Drop-seq barcode designs')
    parser.add_argument('--outdir', help='output dir', required=True)
    parser.add_argument('--adapt', action='append', help='adapter sequence', default=['polyT=A{15}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
    parser.add_argument('--minimum-length', dest='minimum_length', help='minimum_length, default=20', default=20)
    parser.add_argument('--nextseq-trim', dest='nextseq_trim', help='nextseq_trim, default=20', default=20)
    parser.add_argument('--overlap', help='minimum overlap length, default=5', default=5)
    parser.add_argument('--lowQual', type=int, help='max phred of base as lowQual, default=0', default=0)
    parser.add_argument('--lowNum', type=int, help='max number with lowQual allowed, default=2', default=2)

    parser.add_argument('--starMem', help='starMem, default=30', default=30)
    parser.add_argument('--genomeDir', help='genome index dir', required=True)
    parser.add_argument('--refFlat', help='refFlat,for stat mapping region', required=True)
    #parser.add_argument('--runThreadN', type=int, help='', default=2)
    parser.add_argument('--type', help='Specify attribute type in GTF annotation', default='exon')
    parser.add_argument('--annot', help='gtf', required=True)

    parser.add_argument('--cells', type=int, help='cell number, default=3000', default=3000)
    args = vars(parser.parse_args())

    fq_dict, sample_arr, cells_dict = parse_map(args['mapfile'])

    # 链接数据
    raw_dir = args['outdir'] + '/data_give/rawdata'
    os.system('mkdir -p %s'%(raw_dir))
    with open(raw_dir + '/ln.sh', 'w') as fh:
        fh.write('cd %s\n'%(raw_dir))
        for s, arr in fq_dict.items():
            fh.write('ln -sf %s %s\n'%(arr[0], s + '_1.fq.gz'))
            fh.write('ln -sf %s %s\n'%(arr[1], s + '_2.fq.gz'))
    #os.system('sh %s'%(raw_dir+'/ln.sh'))

    logdir = args['outdir']+'/log'
    os.system('mkdir -p %s'%(logdir))
    sjm_cmd = 'log_dir %s\n'%(logdir)
    sjm_order = ''

    for n in sample_arr:
        # barcode
        arr = fq_dict[n]
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='01.barcode')
        cmd = '''source activate scope-tools3; python {app} barcode --fq1 {fq1} --fq2 {fq2} --pattern {pattern} 
                --whitelist {whitelist} --linker {linker} --sample {samplename} --lowQual {lowQual} 
                --lowNum {lowNum} --outdir {outdir};'''.format(
            app = toolsdir + '/scope.py', fq1 = arr[0], fq2 =arr[1], pattern=args['pattern'], 
            whitelist=args['whitelist'], linker=args['linker'], samplename=n, 
            lowQual=args['lowQual'], lowNum=args['lowNum'], outdir=outdir
        )
        sjm_cmd += generate_sjm(cmd, 'barcode_'+n)

        # adapt
        fq = outdir + '/' + n + '_2.fq.gz'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='02.cutadapt')
        cmd = '''source activate scope-tools3; python {app} cutadapt --fq {fq} --sample {samplename} --outdir 
            {outdir}'''.format( app = toolsdir + '/scope.py', fq=fq, samplename = n, outdir = outdir)
        sjm_cmd += generate_sjm(cmd, 'adapt_' + n, m=2)
        sjm_order += 'order adapt_%s after barcode_%s\n'%(n, n)

        # STAR
        fq = outdir + '/' + n + '_clean_2.fq.gz'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='03.STAR')
        cmd = '''source activate scope-tools3; python {app} STAR --fq {fq} --sample {samplename} --refFlat {refFlat} 
        --genomeDir {genomeDir} --thread 8 --outdir {outdir}'''.format(
            app = toolsdir + '/scope.py', fq=fq, samplename=n, refFlat=args['refFlat'], genomeDir=args['genomeDir'],
            outdir = outdir)

        sjm_cmd += generate_sjm(cmd, 'STAR_' + n, m=args['starMem'], x=8)
        sjm_order += 'order STAR_%s after adapt_%s\n'%(n, n)
        
        # featureCounts
        bam = outdir + '/' + n + '_Aligned.sortedByCoord.out.bam'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='04.featureCounts')
        cmd = '''source activate scope-tools3; python {app} featureCounts --input {bam} --annot {annot} --type {type} --sample 
                {samplename} --thread 8 --outdir {outdir}'''.format(
                app = toolsdir + '/scope.py', bam=bam, annot=args['annot'], samplename=n, type=args['type'], outdir = outdir)
        sjm_cmd += generate_sjm(cmd, 'featureCounts_' + n, m=8, x=8)
        sjm_order += 'order featureCounts_%s after STAR_%s\n'%(n, n)

        # count
        bam = outdir + '/' + n + '_name_sorted.bam'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='05.count')
        cmd = '''source activate scope-tools3; python {app} count --bam {bam} --sample {samplename} --cells {cells} 
        --outdir {outdir}'''.format(app=toolsdir + '/scope.py', 
                                            bam=bam, samplename=n, cells =cells_dict[n], outdir=outdir )
        sjm_cmd += generate_sjm(cmd, 'count_' + n, m=30)
        sjm_order += 'order count_%s after featureCounts_%s\n'%(n, n)

    # merged report 
    cmd = '''source activate scope-tools3; python {app} --samples {samples} --workdir {workdir};'''.format(
        app=toolsdir + '/merge_table.py', samples=','.join(sample_arr), workdir=args['outdir'])
    sjm_cmd += generate_sjm(cmd, 'report')
    sjm_order += ''.join(['order report after count_%s\n'%(n) for n in sample_arr])
    

    with open(logdir + '/sjm.job', 'w') as fh:
        fh.write(sjm_cmd+'\n')
        fh.write(sjm_order)

if __name__ == '__main__':
    main()


