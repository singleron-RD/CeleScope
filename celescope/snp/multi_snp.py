from celescope.__init__ import __CONDA__
from celescope.snp.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_snp(Multi):
    def custome_args(self):
        self.STAR_args()
        self.parser.add_argument('--gene_list', help="gene_list", required=True)
        self.parser.add_argument('--probe_file', help="probe fasta file")
        self.parser.add_argument('--annovar_config', help='annovar soft config file')

    def read_custome_args(self):
        self.read_STAR_args()
        self.gene_list = self.args.gene_list
        self.probe_file = self.args.probe_file
        self.annovar_config = self.args.annovar_config
        self.gtf_type = 'gene'
        self.outFilterMatchNmin = 35

    def snpCalling(self, sample):
        step = 'snpCalling'
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_name_sorted.bam'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--bam {bam} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--genomeDir {self.genomeDir} '
            f'--gene_list {self.gene_list} '
            f'--thread {self.thread} '
        )
        self.generate_cmd(cmd, step, sample, m=8, x=self.thread)

    def analysis_snp(self, sample):
        step = 'analysis_snp'
        vcf_anno = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_anno.vcf'
        index_file = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_cell_index.tsv'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--vcf_anno {vcf_anno} '
            f'--index_file {index_file} '
            f'--annovar_config {self.annovar_config} '
        )
        self.generate_cmd(cmd, step, sample, m=8, x=self.thread)


def main():
    multi = Multi_snp(__ASSAY__, __STEPS__, __CONDA__)
    multi.col4_default = None
    multi.run()

if __name__ == '__main__':
    main()

