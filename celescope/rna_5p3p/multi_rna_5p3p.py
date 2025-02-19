from celescope.rna_5p3p.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_rna_5p3p(Multi):
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

    def get_5p3p_fq(self, sample):
        arr = self.fq_dict[sample]
        arr["fq1_3p"] = []
        arr["fq2_3p"] = []
        arr["fq1_5p"] = []
        arr["fq2_5p"] = []
        for fq1, fq2, col4 in zip(arr["fq1"], arr["fq2"], arr["col4"]):
            if col4 not in ("5p", "3p"):
                raise Exception(
                    f"ERROR: The 4th column of mapfile must be 3p or 5p. But {col4} found."
                )
            arr[f"fq1_{col4}"].append(fq1)
            arr[f"fq2_{col4}"].append(fq2)
        return arr

    def convert(self, sample):
        step = "convert"
        arr = self.get_5p3p_fq(sample)
        cmd_line = self.get_cmd_line(step, sample)
        fq1_5p = ",".join(arr["fq1_5p"])
        fq1_3p = ",".join(arr["fq1_3p"])
        fq2_5p = ",".join(arr["fq2_5p"])
        fq2_3p = ",".join(arr["fq2_3p"])
        cmd = (
            f"{cmd_line} "
            f"--fq1_5p {fq1_5p} --fq1_3p {fq1_3p} "
            f"--fq2_5p {fq2_5p} --fq2_3p {fq2_3p} "
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=1)

    def starsolo(self, sample):
        step = "starsolo"
        arr = self.get_5p3p_fq(sample)
        cmd_line = self.get_cmd_line(step, sample)
        fq1_list, fq2_list = [], []
        cnt = {"3p": 0, "5p": 0}
        for x in arr["col4"]:
            cnt[x] += 1
            fq1 = f'{self.outdir_dic[sample]["convert"]}/{sample}_{x}{cnt[x]}_R1.fq.gz'
            fq2 = f'{self.outdir_dic[sample]["convert"]}/{sample}_{x}{cnt[x]}_R2.fq.gz'
            fq1_list.append(fq1)
            fq2_list.append(fq2)
        fq1_str = ",".join(fq1_list)
        fq2_str = ",".join(fq2_list)
        cmd = (
            f"{cmd_line} "
            f"--fq1 {fq1_str} --fq2 {fq2_str} "
            f'--soloFeatures "GeneFull_Ex50pAS Gene SJ" '
        )
        self.process_cmd(cmd, step, sample, m=self.args.starMem, x=self.args.thread)

    def analysis(self, sample):
        step = "analysis"
        matrix_file = f'{self.outdir_dic[sample]["outs"]}/filtered'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--matrix_file {matrix_file} "
        self.process_cmd(cmd, step, sample, m=10, x=1)


def main():
    multi = Multi_rna_5p3p(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
