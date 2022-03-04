from celescope.tcr_fl.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_tcr_fl(Multi):

    def split_fq(self, sample):
        step = 'split_fq'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--nCell {self.nCell} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def assemble(self, sample):
        step = 'assemble'
        fastq_dir = f'{self.outdir_dic[sample]["split_fq"]}/fastq'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fastq_dir {fastq_dir} '
            f'--thread {self.thread} '
        )
        self.process_cmd(cmd, step, sample, m=4 * int(self.args.thread), x=self.args.thread)


def main():
    multi = Multi_tcr_fl(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
