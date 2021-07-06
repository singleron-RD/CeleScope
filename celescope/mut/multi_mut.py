
from celescope.mut.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_mut(Multi):

    def mapping_mut(self, sample):
        step = 'mapping_mut'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--thread {self.thread} '
            f'--indel_genomeDir {self.indel_genomeDir} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def count_mut(self, sample):
        step = 'count_mut'
        bam = f'{self.outdir_dic[sample]["mapping_mut"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step } '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--bam {bam} '
            f'--mut_file {self.mut_file} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--shift_base {self.shift_base} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=1)


def main():
    multi = Multi_mut(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
