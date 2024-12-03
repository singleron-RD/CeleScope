from celescope.pathseq.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_pathseq(Multi):
    """
    ## Usage
    ```
        multi_rna\\
        --mapfile ./rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --thread 8\\
        --mod shell
    ```
    """

    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def pathseq(self, sample):
        step = "pathseq"
        cmd_line = self.get_cmd_line(step, sample)
        input_bam = (
            f'{self.outdir_dic[sample]["outs"]}/{sample}_Aligned.sortedByCoord.out.bam'
        )
        cmd = f"{cmd_line} " f"--input_bam {input_bam} "
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)


def main():
    multi = Multi_pathseq(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
