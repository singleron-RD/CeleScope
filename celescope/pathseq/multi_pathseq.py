from celescope.pathseq.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_pathseq(Multi):
    def starsolo(self, sample):
        step = "starsolo"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f'{cmd_line} ' f'--fq1 {arr["fq1_str"]} --fq2 {arr["fq2_str"]} '
        self.process_cmd(
            cmd,
            step,
            sample,
            m=int(self.args.limitBAMsortRAM / 1e9),
            x=self.args.thread,
        )

    def pathseq(self, sample):
        step = "pathseq"
        cmd_line = self.get_cmd_line(step, sample)
        input_bam = (
            f'{self.outdir_dic[sample]["outs"]}/{sample}_Aligned.sortedByCoord.out.bam'
        )
        cmd = f"{cmd_line} " f"--input_bam {input_bam} "
        self.process_cmd(cmd, step, sample, m=100, x=self.args.thread)

    def count_pathseq(self, sample):
        step = "count_pathseq"
        cmd_line = self.get_cmd_line(step, sample)
        pathseq_bam_file = f"{self.outdir_dic[sample]['pathseq']}/{sample}_pathseq.bam"
        pathseq_score_file = (
            f"{self.outdir_dic[sample]['pathseq']}/{sample}_pathseq_score.txt"
        )
        unmap_bam_file = f"{self.outdir_dic[sample]['pathseq']}/{sample}_unmapped.bam"
        cmd = (
            f"{cmd_line} "
            f"--pathseq_bam_file {pathseq_bam_file} "
            f"--pathseq_score_file {pathseq_score_file} "
            f"--unmap_bam_file {unmap_bam_file} "
            f"--match_dir {self.col4_dict[sample]} "
        )
        self.process_cmd(cmd, step, sample, m=5, x=self.args.thread)

    def analysis_pathseq(self, sample):
        step = "analysis_pathseq"
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f"{cmd_line} "
            f"--umi_matrix_file {self.outdir_dic[sample]['outs']}/{sample}_raw_UMI_matrix.tsv.gz "
            f"--match_dir {self.col4_dict[sample]} "
        )
        self.process_cmd(cmd, step, sample, m=1, x=1)


def main():
    multi = Multi_pathseq(__ASSAY__, min_col=4)
    multi.run()


if __name__ == "__main__":
    main()
