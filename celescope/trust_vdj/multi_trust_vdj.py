from celescope.trust_vdj.__init__ import __ASSAY__
from celescope.tools.Multi import Multi


class Multi_trust_vdj(Multi):

    def trust_assemble(self, sample):
        step = 'trust_assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fq1 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_new_R1.fq{self.fq_suffix}'
        fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_new_R2.fq{self.fq_suffix}' 
        cmd = (
            f'{cmd_line} '
            f'--fq1 {fq1} '
            f'--fq2 {fq2} '
            f'--match_dir {self.col4_dict[sample]}'
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)


    def res_filter(self, sample):
        step = 'res_filter'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_trust_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()
