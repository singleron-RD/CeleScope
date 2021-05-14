from celescope.snp.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_snp(Multi):

    def star(self, sample):
        step = 'star'
        cmd_line = self.get_cmd_line(step, sample)
        if self.args.not_consensus:
            fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        else:
            fq = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fq'
            cmd_line += ' --consensus_fq '
   
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def target_metrics(self, sample):
        step = 'target_metrics'
        cmd_line = self.get_cmd_line(step, sample)
        bam = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam'
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)


    def snpCalling(self, sample):
        step = 'snpCalling'
        cmd_line = self.get_cmd_line(step, sample)
        bam = f'{self.outdir_dic[sample]["target_metrics"]}/{sample}_filtered.bam'
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)

    def analysis_snp(self, sample):
        step = 'analysis_snp'
        vcf = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_merged.vcf'
        CID_file = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_CID.tsv'
        variant_count_file = f'{self.outdir_dic[sample]["snpCalling"]}/{sample}_variant_count.tsv'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--vcf {vcf} '
            f'--CID_file {CID_file} '
            f'--variant_count_file {variant_count_file} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=1)


def main():
    multi = Multi_snp(__ASSAY__)
    multi.run()

if __name__ == '__main__':
    main()

