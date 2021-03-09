from celescope.snp.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_snp(Multi):
    
    def custome_args(self):
        self.STAR_args()
        self.parser.add_argument('--gene_list', help="gene_list", required=True)
        self.parser.add_argument('--annovar_config', help='annovar soft config file')


    def read_custome_args(self):
        self.read_STAR_args()
        self.gene_list = self.args.gene_list
        self.annovar_config = self.args.annovar_config
        self.gtf_type = 'gene'
        self.outFilterMatchNmin = 35

    def STAR(self, sample):
        step = 'STAR'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        if self.args.umi_consensus:
            fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_consensus.fastq.gz'
        if self.outFilterMatchNmin == 0:
            self.outFilterMatchNmin = 35
            
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
    multi.col4_default = None
    multi.run()

if __name__ == '__main__':
    main()

