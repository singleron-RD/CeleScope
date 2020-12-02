from celescope.__init__ import __CONDA__
from celescope.tag.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_tag(Multi):
    def custome_args(self):
        self.parser.add_argument('--thread', help='thread', default=6)
        self.parser.add_argument(
            "--UMI_min",
            help="cells have tag_UMI>=UMI_min are considered as valid cell",
            default="auto")
        self.parser.add_argument("--dim", help="tag dimension", default=1)
        self.parser.add_argument(
            "--SNR_min",
            help="minimum signal to noise ratio",
            default="auto")
        self.parser.add_argument("--fq_pattern", help="tag read2 pattern", required=True)
        self.parser.add_argument("--linker_fasta", help="linker fasta")
        self.parser.add_argument("--barcode_fasta", help="barcode fasta", required=True)

    def read_custome_args(self):
        self.thread = self.args.thread
        self.UMI_min = self.args.UMI_min
        self.dim = self.args.dim
        self.SNR_min = self.args.SNR_min
        self.fq_pattern = self.args.fq_pattern
        self.linker_fasta = self.args.linker_fasta
        self.barcode_fasta = self.args.barcode_fasta

    def mapping_tag(self, sample):
        step = 'mapping_tag'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--fq_pattern {self.fq_pattern} '
            f'--barcode_fasta {self.barcode_fasta} '
            f'--linker_fasta {self.linker_fasta} '
        )
        self.generate_other(cmd, step, sample, m=5, x=1)

    def count_tag(self, sample):
        step = 'count_tag'
        read_count_file = f'{self.outdir_dic[sample]["mapping_tag"]}/{sample}_read_count.tsv'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--read_count_file {read_count_file} '
            f'--dim {self.dim} '
            f'--UMI_min {self.UMI_min} '
            f'--SNR_min {self.SNR_min} '
        )
        self.generate_other(cmd, step, sample, m=5, x=1)


    def analysis_tag(self, sample):
        step = 'analysis_tag'
        tsne_tag_file = f'{self.outdir_dic[sample]["count_tag"]}/{sample}_tsne_tag.tsv'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--tsne_tag_file {tsne_tag_file} '
        )
        self.generate_other(cmd, step, sample, m=5, x=1)



def main():
    multi = Multi_tag(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()

if __name__ == '__main__':
    main()

