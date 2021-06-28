from celescope.trust_vdj.__init__ import __ASSAY__
from celescope.tools.Multi import Multi


class Multi_trust_vdj(Multi):

    def convert(self, sample):
        step = 'convert'
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def matching(self, sample):
        step = 'matching'
        cmd_line = self.get_cmd_line(step, sample)
        fq1 = f'{self.outdir_dic[sample]["convert"]}/{sample}_1.fq{self.fq_suffix}'
        fq2 = f'{self.outdir_dic[sample]["convert"]}/{sample}_2.fq{self.fq_suffix}'
        cmd = (
            f'{cmd_line} '
            f'--fq1 {fq1} '
            f'--fq2 {fq2} '
            f'--match_dir {self.col4_dict[sample]}'
        )
        self.process_cmd(cmd, step, sample, m=5, x=3)


    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fq1 = f'{self.outdir_dic[sample]["matching"]}/{sample}_matched_R1.fq'
        fq2 = f'{self.outdir_dic[sample]["matching"]}/{sample}_matched_R2.fq' 
        cmd = (
            f'{cmd_line} '
            f'--fq1 {fq1} '
            f'--fq2 {fq2} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)


    def mapping(self, sample):
        step = 'mapping'
        cmd_line = self.get_cmd_line(step, sample)
        fq = f'{self.outdir_dic[sample]["assemble"]}//{sample}_toassemble.fq'
        cmd = (
            f'{cmd_line} '
            f'--fq {fq}'
        )
        self.process_cmd(cmd, step, sample, m=5, x=5)


    def res_filter(self, sample):
        step = 'res_filter'
        cmd_line = self.get_cmd_line(step, sample)
        report = f'{self.outdir_dic[sample]["assemble"]}//{sample}_barcode_report.tsv'
        fa = f'{self.outdir_dic[sample]["assemble"]}//{sample}_annot.fa'
        count_file = f'{self.outdir_dic[sample]["matching"]}/count.txt'
        cmd = (
            f'{cmd_line} '
            f'--report {report} '
            f'--fa {fa} '
            f'--count_file {count_file} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_trust_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()
