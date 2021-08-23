from celescope.snp.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_snp(Multi):
    """
    Usage
    ```
    multi_snp\\
        --mapfile ./test1.mapfile\\
        --genomeDir {genomeDir after running celescope snp mkref}\\
        --thread 10\\
        --mod shell\\
        --gene_list gene_list.tsv\\
        --annovar_config annovar.config\\
    ```
    annovar_config file
    ```
    [ANNOVAR]
    dir = /Public/Software/annovar/  
    db = /SGRNJ/Database/script/database/annovar/humandb  
    buildver = hg38  
    protocol = refGene,cosmic70  
    operation = g,f  
    ```
    """

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

    def variant_calling(self, sample):
        step = 'variant_calling'
        cmd_line = self.get_cmd_line(step, sample)
        bam = f'{self.outdir_dic[sample]["target_metrics"]}/{sample}_filtered_sorted.bam'
        cmd = (
            f'{cmd_line} '
            f'--bam {bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)

    def analysis_snp(self, sample):
        step = 'analysis_snp'
        filter_vcf = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_filter.vcf'
        CID_file = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_CID.tsv'
        filter_variant_count_file = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_filter_variant_count.tsv'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--filter_vcf {filter_vcf} '
            f'--CID_file {CID_file} '
            f'--filter_variant_count_file {filter_variant_count_file} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=1)


def main():
    multi = Multi_snp(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
