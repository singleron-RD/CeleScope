from celescope.tracer_vdj.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_tracer_vdj(Multi):
    def custome_args(self):
        self.parser.add_argument('--thread', help='thread', default=20)
        self.parser.add_argument('--mode', help='TCR or BCR', choices=['TCR', 'BCR'])
        self.parser.add_argument('--species', help='species name', choices=['Hsap', 'Mmus'])

    def read_custome_args(self):
        self.thread = self.args.thread
        self.mode = self.args.mode
        self.species = self.args.species

    def split_fastq(self, sample):
        step = 'split_fastq'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--mode {self.mode} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


    def go_assemble(self, sample):
        step = 'go_assemble'
        fastq_dir = f'{self.outdir_dic[sample]["split_fq"]}/fastq'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fastq_dir {fastq_dir} '
            f'--mode {self.mode} '
            f'--species {self.species} '
            f'--thread {self.thread} '
        )
        self.process_cmd(cmd, step, sample, m=1.5 * int(self.args.thread), x=self.args.thread)


def main():
    multi = Multi_tracer_vdj(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()

