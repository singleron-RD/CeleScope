from celescope.tracer_vdj.__init__ import __ASSAY__
from celescope.tools.Multi import Multi


class Multi_tracer_vdj(Multi):

    def split_fastq(self, sample):
        step = 'split_fastq'
        cmd_line = self.get_cmd_line(step, sample)
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
            f'--match_dir {self.col4_dict[sample]}'
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def go_assemble(self, sample):
        step = 'go_assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fastq_dir = f'{self.outdir_dic[sample]["split_fastq"]}/fastq'
        cmd = (
            f'{cmd_line} '
            f'--fastq_dir {fastq_dir} '
        )
        self.process_cmd(cmd, step, sample, m=30, x=self.args.thread)

    def vdj_sum(self, sample):
        step = 'vdj_sum'
        cmd_line = self.get_cmd_line(step, sample)
        ass_dir = f'{self.outdir_dic[sample]["go_assemble"]}'

        fastq_dir = f'{self.outdir_dic[sample]["split_fastq"]}/fastq' 

        cmd = (
            f'{cmd_line} '
            f'--ass_dir {ass_dir} '
            f'--fastq_dir {fastq_dir} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=2)

def main():
    multi = Multi_tracer_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()

