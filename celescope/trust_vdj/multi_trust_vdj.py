from celescope.trust_vdj.__init__ import __ASSAY__
from celescope.tools.multi import Multi

class Multi_trust_vdj(Multi):

    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'

        cmd = (
            f'{cmd_line} '
            f'--fq2 {fq2} '
            f'--match_dir {self.col4_dict[sample]}'
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)
    
    def summarize(self, sample):
        step = 'summarize'
        cmd_line = self.get_cmd_line(step, sample)
        full_len_assembly = f'{self.outdir_dic[sample]["assemble"]}/assemble/{sample}_full_len.fa'
        reads_assignment = f'{self.outdir_dic[sample]["assemble"]}/assemble/{sample}_assign.out'
        assembled_fa = f'{self.outdir_dic[sample]["assemble"]}/assemble/{sample}_assembled_reads.fa'
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--full_len_assembly {full_len_assembly} '
            f'--reads_assignment {reads_assignment} '
            f'--assembled_fa {assembled_fa} '
            f'--fq2 {fq2} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

def main():
    multi = Multi_trust_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()
