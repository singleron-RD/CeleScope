from celescope.tools.multi import Multi
from celescope.trust_vdj.__init__ import __ASSAY__


class Multi_trust_vdj(Multi):

    def convert(self, sample):
        step = 'convert'
        cmd_line = self.get_cmd_line(step, sample)
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--fq2 {fq2} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fq1 = f'{self.outdir_dic[sample]["convert"]}/{sample}_clean_1.fq'
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
        cb_stat = f'{self.outdir_dic[sample]["barcode"]}/stat.txt'
        match_dir = f'{self.col4_dict[sample]}'
        cmd = (
            f'{cmd_line} '
            f'--fq1 {fq1} '
            f'--fq2 {fq2} '
            f'--cb_stat {cb_stat} '
            f'--match_dir {match_dir} '
        )
        self.process_cmd(cmd, step, sample, m=30, x=self.args.thread)
        
    def res_sum(self, sample):
        step = 'res_sum'
        cmd_line = self.get_cmd_line(step, sample)
        all_rep = f'{self.outdir_dic[sample]["assemble"]}/outs/all_contig_annotations.csv'
        fa = f'{self.outdir_dic[sample]["assemble"]}/{sample}_annot.fa'

        cmd = (
            f'{cmd_line} '
            f'--fa {fa} '
            f'--all_rep {all_rep} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_trust_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()
