from celescope.snp.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_snp(Multi):

    
    def snpCalling_args(self):
        self.parser.add_argument('--gene_list', help="gene_list", required=True)

    def read_snpCalling_args(self):
        self.gene_list = self.args.gene_list

    def analysis_snp_args(self):
        self.parser.add_argument('--annovar_config', help='annovar soft config file')

    def read_analysis_snp_args(self):
        self.annovar_config = self.args.annovar_config

    def custome_args(self):
        self.consensus_args()
        self.STAR_args()
        self.snpCalling_args()
        self.analysis_snp_args()
        self.parser.add_argument('--not_consensus', help="not consensus", action="store_true")
        

    def read_custome_args(self):
        self.read_consensus_args()
        self.read_STAR_args()
        self.read_snpCalling_args()
        self.read_analysis_snp_args()
        self.gtf_type = 'gene'
        self.not_consensus = self.args.not_consensus
        self.consensus_fq = not self.not_consensus
        if self.not_consensus:
            self.steps_run = ','.join(['sample', 'barcode', 'cutadapt', 'STAR', 'featureCounts', 'snpCalling', 'analysis_snp'])

    def STAR(self, sample):
        step = 'STAR'
        if self.not_consensus:
            fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        else:
            fq = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fq'         
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--fq {fq} '
            f'--genomeDir {self.genomeDir} '
            f'--thread {self.thread} '
            f'{self.debug_str} '
            f'--outFilterMatchNmin {self.outFilterMatchNmin} '
            f'{self.out_unmapped} '
            f'--STAR_param \"{self.STAR_param}\" '
        )
        self.process_cmd(cmd, step, sample, m=self.starMem, x=self.thread)

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
        self.process_cmd(cmd, step, sample, m=8, x=self.thread)

    def analysis_snp(self, sample):
        step = 'analysis_snp'
        vcf = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_merged.vcf'
        CID_file = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_CID.tsv'
        variant_count_file = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_variant_count.tsv'
        cmd = (
            f'{self.__APP__} '
            f'{self.__ASSAY__} '
            f'{step} '
            f'--outdir {self.outdir_dic[sample][step]} '
            f'--sample {sample} '
            f'--assay {self.__ASSAY__} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--vcf {vcf} '
            f'--CID_file {CID_file} '
            f'--annovar_config {self.annovar_config} '
            f'--variant_count_file {variant_count_file} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=1)


def main():
    multi = Multi_snp(__ASSAY__, __STEPS__)
    multi.run()

if __name__ == '__main__':
    main()

