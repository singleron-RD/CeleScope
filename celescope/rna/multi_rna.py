from celescope.rna.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_rna(Multi):
    """
    ## Usage
    ```
        multi_rna\\
        --mapfile ./rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --thread 8\\
        --mod shell
    ```
    Work for both single cell RNA-Seq and single nuclei RNA-Seq.
    """

    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def analysis(self, sample):
        step = "analysis"
        matrix_file = f'{self.outdir_dic[sample]["outs"]}/filtered'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--matrix_file {matrix_file} "
        self.process_cmd(cmd, step, sample, m=10, x=1)


def main():
    multi = Multi_rna(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
