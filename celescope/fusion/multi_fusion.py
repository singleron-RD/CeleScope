from celescope.__init__ import __CONDA__
from celescope.fusion.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_fusion(Multi):
    def custome_args(self):
        self.STAR_args()
        self.parser.add_argument("--fusion_pos",
        help="first base position of the second gene(0-start),tsv file", required=True)
        self.parser.add_argument("--UMI_min", default=1)

    def read_custome_args(self):
        self.read_STAR_args()
        self.fusion_pos = self.args.fusion_pos
        self.UMI_min = self.args.UMI_min

    def STAR_fusion(self, sample):
        step = 'STAR_fusion'
        input_read = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--input_read {input_read} '
            f'--genomeDir {self.genomeDir} '
            f'--thread {self.thread} '
        )
        self.process_cmd(cmd, step, sample, m=self.starMem, x=self.thread)

    def count_fusion(self, sample):
        step = 'count_fusion'
        bam = f'{self.outdir_dic[sample]["STAR_fusion"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--bam {bam} '
            f'--UMI_min {self.UMI_min} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--fusion_pos {self.fusion_pos} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=1)


def main():
    multi = Multi_fusion(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()


if __name__ == '__main__':
    main()


