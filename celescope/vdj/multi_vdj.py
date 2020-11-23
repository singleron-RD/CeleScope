from celescope.__init__ import __CONDA__
from celescope.vdj.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_vdj(Multi):
    def custome_args(self):
        self.parser.add_argument("--type", help='TCR or BCR', required=True)
        self.parser.add_argument(
        '--iUMI', help='minimum number of UMI of identical receptor type and CDR3', default=1)
        self.parser.add_argument('--thread', help='thread', default=6)

    def read_custome_args(self):
        self.iUMI =  self.args.iUMI
        self.type = self.args.type 
        self.thread = self.args.thread

    def mapping_vdj(self, sample):
        step = 'mapping_vdj'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--type {self.type} '
            f'--thread {self.thread} '
        )
        self.generate_other(cmd, step, sample, m=15, x=self.thread)


    def count_vdj(self, sample):
        # count_vdj
        step = 'count_vdj'
        UMI_count_filter1_file = (
            f'{self.outdir_dic[sample]["mapping_vdj"]}/{sample}_UMI_count_filtered1.tsv'
        )
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--type {self.type} '
            f'--iUMI {self.iUMI} '
            f'--UMI_count_filter1_file {UMI_count_filter1_file} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.generate_other(cmd, step, sample, m=8, x=self.thread)


def main():
    multi = Multi_vdj(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()

if __name__ == '__main__':
    main()

