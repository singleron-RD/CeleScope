from celescope.fusion.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_fusion(Multi):
    """
    Features
    - Generate multi-sample scripts.

    Usage
    ```
    multi_fusion\\
    --mapfile ./fusion.mapfile\\
    --fusion_genomeDir {fusion_genomeDir}\\  
    --mod shell
    ```
    """

    def star_fusion(self, sample):
        step = 'star_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def count_fusion(self, sample):
        step = 'count_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        bam = f'{self.outdir_dic[sample]["star_fusion"]}/{sample}_Aligned.sortedByCoord.out.bam'
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
