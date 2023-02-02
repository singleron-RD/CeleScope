from celescope.fusion.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_fusion(Multi):
    """

    ## Usage
    ```
    multi_fusion\\
    --mapfile ./fusion.mapfile\\
    --fusion_genomeDir {fusion_genomeDir}\\  
    --mod shell\\
    ```

    Use `celescope/data/fusion/{panel}/mkref.sh` to generate the fusion genomeDir. 
    ```
    source activate celescope
    celescope fusion mkref \\
    --genome_name BCR_ABL1_PML_RARA_all \\
    --fasta BCR_ABL1_PML_RARA_all.fasta \\
    --fusion_pos BCR_ABL1_PML_RARA_all_pos.txt\\
    ```

    ## Main Output
    - `05.filter_fusion/{sample}_filtered_UMI.csv`: Filtered fusion UMI counts of each cell barcode.
    """

    def star_fusion(self, sample):
        step = 'star_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def count_fusion(self, sample):
        step = 'count_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        capture_bam = f'{self.outdir_dic[sample]["star_fusion"]}/{sample}_Aligned.sortedByCoord.out.bam'
        cmd = (
            f'{cmd_line} '
            f'--capture_bam {capture_bam} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)

    def filter_fusion(self, sample):
        step = 'filter_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        raw_read_count_file = f'{self.outdir_dic[sample]["count_fusion"]}/{sample}_raw_read_count.json'
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--raw_read_count_file {raw_read_count_file} '
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)

    def analysis_fusion(self, sample):
        step = 'analysis_fusion'
        cmd_line = self.get_cmd_line(step, sample)
        filter_umi_file = f'{self.outdir_dic[sample]["filter_fusion"]}/{sample}_filtered_UMI.csv'
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--filter_umi_file {filter_umi_file} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_fusion(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
