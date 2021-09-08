from celescope.trust_vdj.__init__ import __ASSAY__
from celescope.tools.multi import Multi

class Multi_trust_vdj(Multi):

    # def matching(self, sample):
    #     step = 'matching'
    #     cmd_line = self.get_cmd_line(step, sample)
    #     fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
    #     cmd = (
    #         f'{cmd_line} '
    #         f'--fq2 {fq2} '
    #         f'--match_dir {self.col4_dict[sample]}'
    #     )
    #     self.process_cmd(cmd, step, sample, m=5, x=1)


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


def main():
    multi = Multi_trust_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()
