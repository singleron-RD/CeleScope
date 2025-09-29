from celescope.sweetseq.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_sweetseq(Multi):
    """
    ## Usage
    Before running `multi_sweetseq`, you need to run scRNA-Seq data with CeleScope first.
    ```
    multi_sweetseq \\
        --mapfile ./sweetseq.mapfile\\
        --barcode_fasta celescope/data/sweetseq/sweet_tag_barcode.fasta\\
        --linker_fasta celescope/data/sweetseq/sweet_tag_linker.fasta \\
        --fq_pattern L23C15\\
        --mod shell
    ```
    """

    def mapping_tag(self, sample):
        step = "mapping_tag"
        fq = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--fq {fq} "
        self.process_cmd(cmd, step, sample, m=2, x=1)

    def count_tag(self, sample):
        step = "count_tag"
        read_count_file = (
            f'{self.outdir_dic[sample]["mapping_tag"]}/{sample}_read_count.tsv'
        )
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f"{cmd_line} "
            f"--match_dir {self.col4_dict[sample]} "
            f"--read_count_file {read_count_file} "
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)

    def analysis_tag(self, sample):
        step = "analysis_tag"
        cmd_line = self.get_cmd_line(step, sample)
        tsne_tag_file = f'{self.outdir_dic[sample]["count_tag"]}/{sample}_tsne_tag.tsv'
        cmd = (
            f"{cmd_line} "
            f"--match_dir {self.col4_dict[sample]} "
            f"--tsne_tag_file {tsne_tag_file} "
        )
        self.process_cmd(cmd, step, sample, m=2, x=1)


def main():
    multi = Multi_sweetseq(__ASSAY__, min_col=4)
    multi.run()


if __name__ == "__main__":
    main()
