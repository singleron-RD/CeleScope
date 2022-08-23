from celescope.capture_virus.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_capture_virus(Multi):
    """
    ## Usage
    1. Make a virus genomeDir using `celescope capture_virus mkref`

    2. Generate shell scripts
    ```
    multi_capture_virus \\
    --mapfile {mapfile} \\
    --virus_genomeDir {virus_genomeDir} \\
    --not_consensus \\
    --thread 4 \\
    --mod shell
    ```

    3. Run shell scripts
    ```
    sh ./shell/{sample}.sh
    ```
    """

    def star_virus(self, sample):
        step = 'star_virus'
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

    def count_virus(self, sample):
        step = 'count_virus'
        cmd_line = self.get_cmd_line(step, sample)
        capture_bam = f'{self.outdir_dic[sample]["star_virus"]}/{sample}_virus_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{cmd_line} '
            f'--capture_bam {capture_bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)

    def filter_virus(self, sample):
        step = 'filter_virus'
        cmd_line = self.get_cmd_line(step, sample)
        raw_read_count_file = f'{self.outdir_dic[sample]["count_virus"]}/{sample}_raw_read_count.json'
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--raw_read_count_file {raw_read_count_file} '
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)

    def analysis_virus(self, sample):
        step = 'analysis_virus'
        cmd_line = self.get_cmd_line(step, sample)
        filter_umi_file = f'{self.outdir_dic[sample]["filter_virus"]}/{sample}_filtered_UMI.csv'
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--filter_umi_file {filter_umi_file} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_capture_virus(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
