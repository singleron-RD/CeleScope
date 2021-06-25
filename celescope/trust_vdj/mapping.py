import pandas as pd
import glob
from celescope.tools.Step import Step, s_common
from celescope.tools import utils
import os
import re


class Mapping(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        self.outdir = args.outdir
        self.match_dir = args.match_dir
        self.Seqtype = args.Seqtype
        self.sample = args.sample
        self.species = args.species

    @utils.add_log
    def align(self):
        species = self.species
        outdir = self.outdir
        Seqtype = self.Seqtype
        
        stat_file = self.outdir + '/stat.txt'
        fq = f'{outdir}/../02.assemble/TRUST4/{self.sample}_toassemble.fq'

        mapping_summary = []

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
                        mapping_summary.append({
                            'item': item,
                            'count': count,
                            'total_count': total_count,
                        })
        # os.system(f'rm {outdir}/{locus}.sam')   

        # total mapping
        cmd = (
                f'source activate full_len_VDJ; '
                f'bowtie2 -p 5 -k 1 --np 0 --rdg 1,1 --rfg 1,1 '
                f'-x /SGRNJ03/randd/zhouxin/software/TRUST4/index/{species}/{Seqtype} '
                f'-U {fq} '
                f'-S {outdir}/{Seqtype}.sam > {outdir}/log 2>&1'        
        )
        os.system(cmd)
        with open(f'{outdir}/log') as fh: 
            for line in fh:
                if 'reads; of these:' in line:
                    attr = re.findall(r'\d+', line)
                    total_count = int(attr[0])
                if 'aligned exactly 1 time' in line:
                    res = re.findall(r"\d+", line)
                    count = int(res[0])
                    mapping_summary.insert(0, {
                        'item': stat_string,
                        'count': count,
                        'total_count': total_count,
                    })

        os.system(f'rm {outdir}/*.sam')
        os.system(f'rm {outdir}/log')

        df = pd.DataFrame(mapping_summary, columns=['item', 'count', 'total_count'])

        utils.gen_stat(df, stat_file)

    @utils.add_log
    def run(self):
        self.align()

        self.clean_up()

    
def mapping(args):
    step_name = 'mapping'
    mapping_obj = Mapping(args, step_name)
    mapping_obj.run()


def get_opts_mapping(parser, sub_program):
    if sub_program:
        parser = s_common(parser)

    parser.add_argument('--Seqtype', help='select TCR or BCR', choices=["TCR", "BCR"], required=True)
    parser.add_argument('--species', help='species', choices=["Mmus", "Hsap"], required=True)


