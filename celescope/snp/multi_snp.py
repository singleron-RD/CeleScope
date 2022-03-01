from celescope.snp.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_snp(Multi):
    """
    ## Usage

    ### Make a snp reference genomeDir

    1. Run `celescope rna mkref`. If you already have a rna genomeDir, you can use it and skip this step.
    2. Run `celescope snp mkref` under the rna genomeDir. Check [mkref.md](./mkref.md) for help.

    ### Install ANNOVAR, download the annotation database and write a annovar config file.
    https://annovar.openbioinformatics.org/en/latest/user-guide/download/

    ```
    perl /Public/Software/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar cosmic70 humandb/
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

    ### Run multi_snp
    There are two ways to run `multi_snp`

    1. Do not perform consensus before alignment and report read count(recommended for data generated with FocuSCOPE kit).

    ```
    multi_snp\\
        --mapfile ./test1.mapfile\\
        --genomeDir {genomeDir after running celescope snp mkref}\\
        --thread 4\\
        --mod shell\\
        --panel lung_1\\
        --annovar_config annovar.config\\
        --not_consensus
    ```

    2. Do consensus before alignment and report UMI count. 

    ```
    multi_snp\\
        --mapfile ./test1.mapfile\\
        --genomeDir {genomeDir after running celescope snp mkref}\\
        --thread 4\\
        --mod shell\\
        --panel lung_1\\
        --annovar_config annovar.config\\
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
            f'--add_RG '
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
        self.process_cmd(cmd, step, sample, m=8, x=1)

    def filter_snp(self, sample):
        step ='filter_snp'
        vcf = f'{self.outdir_dic[sample]["variant_calling"]}/{sample}_norm.vcf'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--vcf {vcf} '
        )
        self.process_cmd(cmd, step, sample, m=1, x=1)

    def analysis_snp(self, sample):
        step = 'analysis_snp'
        vcf = f'{self.outdir_dic[sample]["filter_snp"]}/{sample}_filtered.vcf'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--vcf {vcf} '
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)


def main():
    multi = Multi_snp(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
