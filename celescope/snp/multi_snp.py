from celescope.snp.__init__ import __STEPS__, __ASSAY__
from celescope.tools.Multi import Multi


class Multi_snp(Multi):

    def STAR(self, sample):
        step = 'STAR'
        cmd_line = self.get_cmd_line(step, sample)
        if self.args.not_consensus:
            fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq.gz'
        else:
            fq = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fq'         
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

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
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)

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
    multi = Multi_snp(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()

