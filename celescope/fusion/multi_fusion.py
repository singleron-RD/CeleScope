from celescope.fusion.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_fusion(Multi):

    def STAR_fusion(self, sample):
        step = 'STAR_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        input_read = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{cmd_line} '
            f'--input_read {input_read} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def count_fusion(self, sample):
        step = 'count_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        bam = f'{self.outdir_dic[sample]["STAR_fusion"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=1)


def main():
    multi = Multi_fusion(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()


