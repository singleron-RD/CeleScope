from celescope.rna.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_rna(Multi):
    """
    Usage
    ```
        multi_rna\\
        --mapfile ./rna.mapfile\\
        --genomeDir /SGRNJ/Public/Database/genome/homo_mus\\
        --thread 8\\
        --mod shell
    ```

    If Single nuclei RNA-Seq is used, you need to add `--gtf_type gene` to include reads mapped to 
    intronic regions.
    """


def main():
    multi = Multi_rna(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
