from celescope.bulk_dynaseq.__init__ import __ASSAY__
from celescope.bulk_rna.multi_bulk_rna import Multi_bulk_rna
from celescope.tools.__init__ import (
    FILTERED_MATRIX_DIR_SUFFIX,
    BARCODE_FILE_NAME,
    OUTS_DIR,
)


class Multi_bulk_dynaseq(Multi_bulk_rna):
    """
    ## Usage

    ```
        multi_bulk_dynaseq\\
        --mapfile ./bulk_dynaseq.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --well_sample ./well_sample.tsv
    ```

    For control sample, set --control to skip replacement step.
    ```
        multi_bulk_dynaseq\\
        --mapfile ./bulk_dynaseq.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --well_sample ./well_sample.tsv\\
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
            f'--SAM_attributes MD '
            f'--STAR_param "--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3" '
        )
        self.process_cmd(
            cmd,
            step,
            sample,
            m=int(self.args.limitBAMsortRAM / 1e9),
            x=self.args.thread,
        )

    def split(self, sample):
        step = "split"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr["fq1_str"]} '
            f'--bam {self.outdir_dic[sample][OUTS_DIR]}/{sample}_Aligned.sortedByCoord.out.bam '
            f'--barcodes {self.outdir_dic[sample][OUTS_DIR]}/filtered/barcodes.tsv.gz '
            f'--split_bam '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def conversion(self, sample):
        step = "conversion"
        split_dir = self.outdir_dic[sample]["split"]
        barcodes = f'{self.outdir_dic[sample]["outs"]}/{FILTERED_MATRIX_DIR_SUFFIX}/{BARCODE_FILE_NAME}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--bam_dir {split_dir} " f"--barcodes {barcodes} "
        self.process_cmd(
            cmd, step, sample, m=self.args.conversionMem, x=self.args.thread
        )

    def substitution(self, sample):
        step = "substitution"
        conversion_dir = self.outdir_dic[sample]["conversion"]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--bam_dir {conversion_dir}"
        self.process_cmd(cmd, step, sample, m=1, x=1)

    def replacement(self, sample):
        step = "replacement"
        conversion_dir = self.outdir_dic[sample]["conversion"]
        cmd_line = self.get_cmd_line(step, sample)
        bg_para = ""
        if sample in self.col5_dict:
            bg_para = f"--bg {self.col5_dict[sample]} "
        cmd = f"{cmd_line} " f"--bam_dir {conversion_dir} " f"{bg_para} "
        self.process_cmd(
            cmd, step, sample, m=self.args.replacementMem, x=self.args.thread
        )


def main():
    multi = Multi_bulk_dynaseq(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
