from celescope.__init__ import __CONDA__
from celescope.capture_virus.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_capture_virus(Multi):
    def custome_args(self):
        self.STAR_args()
        self.parser.add_argument('--virus_genomeDir', help='virus_genomeDir', required=True)

    def read_custome_args(self):
        self.read_STAR_args()
        self.virus_genomeDir = self.args.virus_genomeDir

    def STAR_virus(self, sample):
        step = 'STAR_virus'
        input_read = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--input_read {input_read} '
            f'--virus_genomeDir {self.virus_genomeDir} '
            f'--thread {self.thread} '
        )
        self.process_cmd(cmd, step, sample, m=self.starMem, x=self.thread)

    def count_capture_virus(self, sample):
        step = 'count_capture_virus'
        virus_bam = f'{self.outdir_dic[sample]["STAR_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--virus_bam {virus_bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def analysis_capture_virus(self, sample):        
        step = 'analysis_capture_virus'
        virus_file = f'{self.outdir_dic[sample]["count_capture_virus"]}/{sample}_virus_UMI_count.tsv'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--virus_file {virus_file} '
            f'--match_dir {self.col4_dict[sample]} '

        )
        self.process_cmd(cmd, step, sample, m=15, x=1)


def main():
    multi = Multi_capture_virus(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()


if __name__ == '__main__':
    main()



