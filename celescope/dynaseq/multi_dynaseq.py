from celescope.dynaseq.__init__ import __ASSAY__
from celescope.rna.multi_rna import Multi_rna
from celescope.tools.__init__ import (
    FILTERED_MATRIX_DIR_SUFFIX,
    BARCODE_FILE_NAME,
    FEATURE_FILE_NAME,
    STARSOLO_BAM_SUFFIX,
    OUTS_DIR,
)


class Multi_dynaseq(Multi_rna):
    """
    ## Usage

    ```
        multi_dynaseq\\
        --mapfile ./rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
    ```

    For control sample, set --control to skip replacement step.
    ```
        multi_dynaseq\\
        --mapfile ./rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --control
    ```
    """

    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
            f'--STAR_param "--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3" '
            f'--SAM_attributes MD '
        )
        self.process_cmd(
            cmd,
            step,
            sample,
            m=int(self.args.limitBAMsortRAM / 1e9),
            x=self.args.thread,
        )

    def conversion(self, sample):
        step = "conversion"
        bam = f'{self.outdir_dic[sample]["outs"]}/{sample}_{STARSOLO_BAM_SUFFIX}'
        star_log = f'{self.outdir_dic[sample]["starsolo"]}/{sample}_Log.final.out'
        cell = f'{self.outdir_dic[sample]["outs"]}/{FILTERED_MATRIX_DIR_SUFFIX}/{BARCODE_FILE_NAME}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--bam {bam} " f"--cell {cell} " f"--star_log {star_log} "
        self.process_cmd(
            cmd, step, sample, m=self.args.conversionMem, x=self.args.thread
        )

    def substitution(self, sample):
        step = "substitution"
        bam = f'{self.outdir_dic[sample]["conversion"]}/{sample}.PosTag.bam'
        snp = f'{self.outdir_dic[sample]["conversion"]}/{sample}.snp.csv'
        cmd_line = self.get_cmd_line(step, sample)
        bg_para = ""
        if sample in self.col5_dict:
            bg_para = f"--bg {self.col5_dict[sample]} "
        cmd = f"{cmd_line} " f"--bam {bam} " f"--bg {snp} {bg_para} "
        self.process_cmd(cmd, step, sample, m=1, x=1)

    def replacement(self, sample):
        step = "replacement"
        bam = f'{self.outdir_dic[sample]["conversion"]}/{sample}.PosTag.bam'
        snp = f'{self.outdir_dic[sample]["conversion"]}/{sample}.snp.csv'
        tsne_file = f"{self.outdir_dic[sample][OUTS_DIR]}/tsne_coord.tsv"
        cell = f'{self.outdir_dic[sample]["outs"]}/{FILTERED_MATRIX_DIR_SUFFIX}/{BARCODE_FILE_NAME}'
        gene = f'{self.outdir_dic[sample]["outs"]}/{FILTERED_MATRIX_DIR_SUFFIX}/{FEATURE_FILE_NAME}'
        cmd_line = self.get_cmd_line(step, sample)
        bg_para = ""
        if sample in self.col5_dict:
            bg_para = f"--bg {self.col5_dict[sample]} "
        cmd = (
            f"{cmd_line} "
            f"--bam {bam} "
            f"--bg {snp} {bg_para} "
            f"--tsne {tsne_file} "
            f"--cell {cell} "
            f"--gene {gene} "
        )
        self.process_cmd(
            cmd, step, sample, m=self.args.replacementMem, x=self.args.thread
        )


def main():
    multi = Multi_dynaseq(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
