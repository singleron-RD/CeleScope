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


def main():
    multi = Multi_rna(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
