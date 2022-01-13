from celescope.tag.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_tag(Multi):
    """
    `multi_tag` is used to analyze Single-cell Multiplexing data generated with CLindex<sup>TM</sup> Sample Multiplexing kits. 
    Before running `multi_tag`, you need to run scRNA-Seq data with CeleScope first.

    Usage

    ```
    multi_tag \\
        --mapfile ./tag.mapfile\\
        --barcode_fasta ./tag_barcode.fasta\\
        --fq_pattern L25C15\\
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

    def count_tag(self, sample):
        step = 'count_tag'
        read_count_file = f'{self.outdir_dic[sample]["mapping_tag"]}/{sample}_read_count.tsv'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--read_count_file {read_count_file} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def analysis_tag(self, sample):
        step = 'analysis_tag'
        tsne_tag_file = f'{self.outdir_dic[sample]["count_tag"]}/{sample}_tsne_tag.tsv'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--tsne_tag_file {tsne_tag_file} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def split_tag(self, sample):
        step = 'split_tag'
        umi_tag_file = f'{self.outdir_dic[sample]["count_tag"]}/{sample}_umi_tag.tsv'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--umi_tag_file {umi_tag_file} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_tag(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
