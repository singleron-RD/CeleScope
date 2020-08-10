#!/bin/env python
#coding=utf8

import os
import glob
import sys
import argparse
import re
from collections import defaultdict

toolsdir = os.path.realpath(sys.path[0] + '/../../tools')
multidir = os.path.realpath(sys.path[0] + '/../multi_samples')

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


def parse_map(mapfile, cells):
    fq_dict = defaultdict(list)
    cells_dict = defaultdict(list)
    with open(mapfile) as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue
            tmp = line.split()
            library_id = tmp[0]
            library_path = tmp[1]
            sample_name = tmp[2]
            if len(tmp) == 4:
                cell_number = int(tmp[3])
            else:
                cell_number = cells
            try:
                pattern1_1 = library_path + '/' + library_id + '*' + '_1.fq.gz'
                pattern1_2 = library_path + '/' + library_id + '*' + 'R1_*.fastq.gz'
                pattern2_1 = library_path + '/' + library_id + '*' + '_2.fq.gz'
                pattern2_2 = library_path + '/' + library_id + '*' + 'R2_*.fastq.gz'
                fq1 = (glob.glob(pattern1_1) + glob.glob(pattern1_2))[0]
                fq2 = (glob.glob(pattern2_1) + glob.glob(pattern2_2))[0]
            except IndexError as e:
                sys.exit("Error:"+str(e))
                
            assert os.path.exists(fq1), '%s not exists!'%(fq1)
            assert os.path.exists(fq2), '%s not exists!'%(fq2)
            if sample_name in fq_dict:
                fq_dict[sample_name][0].append(fq1)
                fq_dict[sample_name][1].append(fq2)
            else:
                fq_dict[sample_name] = [[fq1], [fq2]]
            cells_dict[sample_name] = cell_number
    
    for sample_name in fq_dict:
        fq_dict[sample_name][0] = ",".join(fq_dict[sample_name][0])
        fq_dict[sample_name][1] = ",".join(fq_dict[sample_name][1])

    return fq_dict, cells_dict


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
    parser.add_argument('--pattern', help='read1 pattern', default='C8L16C8L16C8U8T18')
    parser.add_argument('--bcType', help='choice of barcode types. Currently support scope and Drop-seq barcode designs')
    parser.add_argument('--outdir', help='output dir', default="./")
    parser.add_argument('--adapt', action='append', help='adapter sequence', default=['polyT=A{15}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'])
    parser.add_argument('--minimum-length', dest='minimum_length', help='minimum_length', default=20)
    parser.add_argument('--nextseq-trim', dest='nextseq_trim', help='nextseq_trim', default=20)
    parser.add_argument('--overlap', help='minimum overlap length, default=5', default=5)
    parser.add_argument('--lowQual', type=int, help='max phred of base as lowQual', default=0)
    parser.add_argument('--lowNum', type=int, help='max number with lowQual allowed', default=2)
    parser.add_argument('--starMem', help='starMem', default=30)
    parser.add_argument('--genomeDir', help='genome index dir', required=True)
    parser.add_argument('--gtf_type', help='Specify attribute type in GTF annotation, default=exon', default='exon')
    parser.add_argument('--cells', type=int, help='cell number, default=3000', default=3000)
    parser.add_argument('--conda', help='conda env name', default="scope1.0")
    args = vars(parser.parse_args())

    fq_dict, cells_dict = parse_map(args['mapfile'],args['cells'])

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
    conda = args['conda']

    for n in fq_dict:
        # sample
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='00.sample')
        cmd = '''source activate {conda}; python {app} sample 
        --sample {samplename} --outdir {outdir} --genomeDir {genomeDir};'''.format(
            conda=conda, app=toolsdir + '/scope.py', samplename=n, outdir=outdir, genomeDir=args['genomeDir'])
        sjm_cmd += generate_sjm(cmd, 'sample_'+n)

        # barcode
        arr = fq_dict[n]
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='01.barcode')
        cmd = '''source activate {conda}; python {app} barcode --fq1 {fq1} --fq2 {fq2} --pattern {pattern} 
                --whitelist {whitelist} --linker {linker} --sample {samplename} --lowQual {lowQual} 
                --lowNum {lowNum} --outdir {outdir} --thread 2;'''.format(
            conda=conda, app=toolsdir + '/scope.py', fq1 = arr[0], fq2 =arr[1], pattern=args['pattern'], 
            whitelist=args['whitelist'], linker=args['linker'], samplename=n, 
            lowQual=args['lowQual'], lowNum=args['lowNum'], outdir=outdir
        )
        sjm_cmd += generate_sjm(cmd, 'barcode_'+n, m=5, x=2)
        sjm_order += 'order barcode_%s after sample_%s\n'%(n, n)

        # adapt
        fq = outdir + '/' + n + '_2.fq.gz'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir = args['outdir'], sampledir = n, step='02.cutadapt')
        cmd = '''source activate {conda}; python {app} cutadapt --fq {fq} --sample {samplename} --outdir 
            {outdir}'''.format(conda=conda, app=toolsdir + '/scope.py', fq=fq, samplename=n, outdir=outdir)
        sjm_cmd += generate_sjm(cmd, 'adapt_' + n, m=2)
        sjm_order += 'order adapt_%s after barcode_%s\n'%(n, n)

        # STAR
        fq = outdir + '/' + n + '_clean_2.fq.gz'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir=args['outdir'], sampledir=n, step='03.STAR')
        cmd = '''source activate {conda}; python {app} STAR --fq {fq} --sample {samplename} 
        --genomeDir {genomeDir} --thread 8 --outdir {outdir}'''.format(
            conda=conda, app=toolsdir + '/scope.py', fq=fq, samplename=n, genomeDir=args['genomeDir'],
            outdir=outdir)

        sjm_cmd += generate_sjm(cmd, 'STAR_' + n, m=args['starMem'], x=8)
        sjm_order += 'order STAR_%s after adapt_%s\n'%(n, n)
        
        # featureCounts
        bam = outdir + '/' + n + '_Aligned.sortedByCoord.out.bam'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir=args['outdir'], sampledir=n, step='04.featureCounts')
        cmd = '''source activate {conda}; python {app} featureCounts --input {bam} --type {gtf_type} --sample 
                {samplename} --thread 8 --outdir {outdir}'''.format(
                conda=conda, app=toolsdir + '/scope.py', bam=bam,
                samplename=n, gtf_type=args['gtf_type'], outdir=outdir)
        sjm_cmd += generate_sjm(cmd, 'featureCounts_' + n, m=8, x=8)
        sjm_order += 'order featureCounts_%s after STAR_%s\n'%(n, n)

        # count
        bam = outdir + '/' + n + '_name_sorted.bam'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir=args['outdir'], sampledir=n, step='05.count')
        cmd = '''source activate {conda}; python {app} count --bam {bam} --sample {samplename} --cells {cells} 
            --outdir {outdir}'''.format(conda=conda, app=toolsdir + '/scope.py',  
            bam=bam, samplename=n, cells=cells_dict[n], outdir=outdir)
        sjm_cmd += generate_sjm(cmd, 'count_' + n, m=30)
        sjm_order += 'order count_%s after featureCounts_%s\n'%(n, n)

        # analysis
        matrix_file = outdir + '/' + n + '_matrix.xls'
        outdir = '{basedir}/{sampledir}/{step}'.format(basedir=args['outdir'], sampledir=n, step='06.analysis')
        cmd = '''source activate {conda}; python {app} analysis --matrix_file {matrix_file} --sample {samplename}  
            --outdir {outdir} '''.format(conda=conda, app=toolsdir + '/scope.py', 
            matrix_file=matrix_file, samplename=n,  outdir=outdir)
        sjm_cmd += generate_sjm(cmd, 'analysis_' + n, m=10)
        sjm_order += 'order analysis_%s after count_%s\n'%(n, n)


    # merged report 
    cmd = '''source activate {conda}; python {app} --samples {samples} --workdir {workdir};'''.format(
        conda=conda, app=multidir + '/merge_table.py', samples=','.join(fq_dict.keys()), workdir=args['outdir'])
    sjm_cmd += generate_sjm(cmd, 'report')
    sjm_order += ''.join(['order report after count_%s\n'%(n) for n in fq_dict])
    with open(logdir + '/sjm.job', 'w') as fh:
        fh.write(sjm_cmd+'\n')
        fh.write(sjm_order)


if __name__ == '__main__':
    main()


