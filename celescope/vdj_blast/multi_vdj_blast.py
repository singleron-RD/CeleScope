from celescope.tools.multi import Multi
from celescope.vdj_blast.__init__ import __ASSAY__

class Multi_vdj_blast(Multi):
    """
    ## Usage
    ```
    multi_vdj_blast \\
        --mapfile ./vdj.mapfile \\
        --thread 8 \\
        --mod shell
    ``` 
    """
    def mapping_vdj(self, sample):
        step = 'mapping_vdj'
        cmd_line = self.get_cmd_line(step, sample)
        fasta = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fasta'
        cmd = (
            f'{cmd_line} '
            f'--fasta {fasta} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)


def main():
    multi = Multi_vdj_blast(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()