from celescope.bulk_rna.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_bulk_rna(Multi):
    """
    ## Usage
    - For default 384 wells:
    ```
        multi_bulk_rna\\
        --mapfile ./bulk_rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus
    ```
    - For 96 wells, set `--wells`:
    ```
        multi_bulk_rna\\
        --mapfile ./bulk_rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --wells 96
    ```

    Work for AccuraCode bulk RNA-Seq.
    """

    def count(self, sample):
        step = 'count'
        count_detail = f'{self.outdir_dic[sample]["featureCounts"]}/{sample}_count_detail.txt'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--count_detail {count_detail} '
        )
        self.process_cmd(cmd, step, sample, m=20, x=1)


def main():
    multi = Multi_bulk_rna(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
