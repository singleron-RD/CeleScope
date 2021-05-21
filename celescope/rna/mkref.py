import subprocess

import celescope.tools.utils as utils
from celescope.tools.step import Step, s_common


class Mkref(Step):
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        self.genome_name = self.args.genome_name
        self.fasta = self.args.fasta
        self.gtf = self.args.gtf
        self.mt_gene_list = self.args.mt_gene_list

    @utils.add_log
    def build_star_index(self):
        cmd = f'''STAR \
    --runMode genomeGenerate \
    --runThreadN {self.thread} \
    --genomeDir {self.outdir} \
    --genomeFastaFiles {self.fasta} \
    --sjdbGTFfile {self.gtf} \
    --sjdbOverhang 100'''
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    
    def run(self):
        self.build_star_index()

def mkref(args):
    step_name = 'mkref'
    runner = Mkref(args, step_name)
    runner.run()

def get_opts_mkref(self):
    

    

        
