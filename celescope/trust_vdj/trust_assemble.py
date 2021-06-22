import os
from celescope.tools import utils
from celescope.tools.Step import Step, s_common
from celescope.tracer_vdj.split_fastq import get_barcodes
from celescope.tools.barcode import *
import pysam
import pandas as pd
from collections import defaultdict


TRUST = '/SGRNJ03/randd/zhouxin/software/TRUST4/run-trust4'


def count_fq(fq1):
    dic = defaultdict(list)
    with pysam.FastxFile(fq1) as fq:
        for entry in fq:
            attr = entry.sequence
            cb = attr[:24]
            umi = attr[24:]
            name = entry.name
            dic['barcode'].append(cb)
            dic['UMI'].append(umi)
            dic['seq_name'].append(name)

    count_df = pd.DataFrame(dic, columns=list(dic.keys()))
    
    return count_df


@utils.add_log
def match_barcodes(outdir, match_dir, Seqtype, fq1):
    annotated_bcs = get_barcodes(match_dir, Seqtype)
    bcs_df = pd.DataFrame(annotated_bcs, columns=['barcode'])
    count_df = count_fq(fq1)

    # count UMI
    df_umi = count_df.groupby(['barcode', 'UMI'], as_index=False).agg({'seq_name': 'count'})
    df_umi = df_umi.groupby(['barcode'], as_index=False).agg({'UMI': 'count'})
    df_umi = df_umi.sort_values(by='UMI', ascending=False)
    df_umi.to_csv(f'{outdir}/count.txt', sep='\t', index=False)

    df_n = pd.merge(bcs_df, count_df, on='barcode', how='inner')
    seqnames = df_n['seq_name'].tolist()
    seqlist = open(f'{outdir}/seqlist.txt', 'w')
    for name in seqnames:
        seqlist.write(str(name) + '\n')


def mapping_summary(outdir, Seqtype, fq, species):
    
    stat_file = outdir + '/stat.txt'

    trust_assemble_summary = []

    total_mapped = 0

    #with pysam.FastxFile(fq) as fh:
        #total_count = 0
        #for entry in fh:
            #total_count += 1

    if Seqtype == 'TCR':
        loci = ['TRA', 'TRB']
        stat_string = 'All reads Mapped to TRA and TRB' 

    elif Seqtype == 'BCR':
        loci = ['IGH', 'IGL', 'IGK']
        stat_string = 'All reads Mapped to IGH, IGL and IGK'

    for locus in loci:
        cmd = (
            f'source activate bracer; '
            f'bowtie2 -p 5 -k 1 --np 0 --rdg 1,1 --rfg 1,1 '
            f'-x /SGRNJ03/randd/zhouxin/software/TRUST4/index/{species}/{locus} '
            f'-U {fq} '
            f'-S {outdir}/{locus}.sam > {outdir}/log 2>&1'
        )
        os.system(cmd)

        with open(f'{outdir}/log') as fh:
            for line in fh:
                if 'reads; of these:' in line:
                    attr = re.findall(r'\d+', line)
                    total_count = int(attr[0])
                if 'aligned exactly 1 time' in line:
                    res = re.findall(r"\d+", line)
                    item = f'Reads mapped to {locus}'
                    count = int(res[0])
                    total_mapped += count
                    trust_assemble_summary.append({
                        'item': item,
                        'count': count,
                        'total_count': total_count,
                    })

        os.system(f'rm {outdir}/{locus}.sam')

    trust_assemble_summary.insert(0, {
        'item': stat_string,
        'count': total_mapped,
        'total_count': total_count
    })

    os.system(f'rm {outdir}/log')

    df = pd.DataFrame(trust_assemble_summary, columns=['item', 'count', 'total_count'])

    utils.gen_stat(df, stat_file)


class Trust_assemble(Step):
    """
    Features

    - Get fq file
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.match_dir = args.match_dir
        self.Seqtype = args.Seqtype
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.sample = args.sample
        self.species = args.species
        self.speed_up = args.speed_up

    
    @utils.add_log
    def getFqfile(self):
        match_barcodes(self.outdir, self.match_dir, self.Seqtype, self.fq1)

        cmd1 = (
            f'seqtk subseq {self.fq1} {self.outdir}/seqlist.txt > {self.outdir}/{self.sample}_matched_R1.fq'
        )
        os.system(cmd1)

        cmd2 = (
            f'seqtk subseq {self.fq2} {self.outdir}/seqlist.txt > {self.outdir}/{self.sample}_matched_R2.fq'
        )
        os.system(cmd2)


    @utils.add_log
    def run(self):

        self.getFqfile()

        species = self.species

        index_file = f'/SGRNJ03/randd/zhouxin/software/TRUST4/index/{species}/{species}_ref.fa'
        ref = f'/SGRNJ03/randd/zhouxin/software/TRUST4/index/{species}/{species}_IMGT+C.fa'

        string1 = ''
        if self.speed_up:
            string1 = '--repseq '
        cmd = (
            f'{TRUST} -t {self.thread} '
            f'-u {self.outdir}/{self.sample}_matched_R2.fq '
            f'--barcode {self.outdir}/{self.sample}_matched_R1.fq '
            f'--barcodeRange 0 23 + '
            f'-f {index_file} '
            f'--ref {ref} '
            f'{string1}'
            f'-o {self.sample} --od {self.outdir}/TRUST4' 
        )

        Trust_assemble.run.logger.info(cmd)

        if not os.path.exists(f'{self.outdir}/TRUST4/{self.sample}_barcode_report.tsv'):
            os.system(cmd)

            #fq = f'{self.outdir}/TRUST4/{self.sample}_toassemble.fq'

        mapping_summary(self.outdir, self.Seqtype, self.fq2, species)

        os.remove(f'{self.outdir}/seqlist.txt')

        self.clean_up()


@utils.add_log
def trust_assemble(args):
    step_name = 'trust_assemble'
    trust_assemble_obj = Trust_assemble(args, step_name)
    trust_assemble_obj.run()


def get_opts_trust_assemble(parser, sub_program):
    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--fq1', help='R1 reads from barcode step', required=True)
        parser.add_argument('--fq2', help='R2 reads from barcode step', required=True)
        parser.add_argument('--match_dir', help='match_dir', required=True)
    parser.add_argument('--Seqtype', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)
    parser.add_argument('--speed_up', help='speed assemble for TCR/BCR seq data', action='store_true')       








