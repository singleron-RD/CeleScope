
from celescope.citeseq.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_citeseq(Multi):
    """
    ## Usage
    ```
    multi_citeseq \\
        --mapfile ./test.mapfile \\
        --barcode_fasta ./CLindex_TAG.fasta \\
        --fq_pattern L25C15 \\
        --mod shell
    ```
    """

    def mapping_tag(self, sample):
        step = 'mapping_tag'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def count_cite(self, sample):

        step = 'count_cite'
        cmd_line = self.get_cmd_line(step, sample)
        read_count_file = f'{self.outdir_dic[sample]["mapping_tag"]}/{sample}_read_count.tsv'
        cmd = (
            f'{cmd_line} '
            f'--read_count_file {read_count_file} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=1, x=1)

    def analysis_cite(self, sample):

        step = 'analysis_cite'
        cmd_line = self.get_cmd_line(step, sample)
        tsne_coord = f'{self.outdir_dic[sample]["count_cite"]}/{sample}_filtered_tsne_coord.tsv'
        cmd = (
            f'{cmd_line} '
            f'--tsne_coord {tsne_coord} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_citeseq(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
