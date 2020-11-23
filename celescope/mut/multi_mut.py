from celescope.__init__ import __CONDA__
from celescope.mut.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_mut(Multi):
    
    def custome_args(self):
        self.STAR_args()
        self.parser.add_argument("--mut_file", help="mutation file", required=True)
        self.parser.add_argument("--shift_base", default=2)
        self.parser.add_argument('--indel_genomeDir', help='insertion or deletion STAR indexed genome directory', required=True)

    def read_custome_args(self):
        self.read_STAR_args()
        self.mut_file = self.args.mut_file
        self.shift_base = self.args.shift_base
        self.indel_genomeDir = self.args.indel_genomeDir    

    def mapping_mut(self, sample):
        step = 'mapping_mut'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
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
        self.generate_other(cmd, step, sample, m=self.starMem, x=self.thread)

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
            f'shift_base {self.shift_base} '
        )
        self.generate_other(cmd, step, sample, m=8, x=1)


def main():
    multi = Multi_mut(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()

if __name__ == '__main__':
    main()

