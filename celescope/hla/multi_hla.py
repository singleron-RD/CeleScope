from celescope.tools.multi import Multi
from celescope.fusion.__init__ import __ASSAY__


class Multi_hla(Multi):
    """
    Usage
    ```
    multi_hla\\
    --mapfile ./hla.mapfile\\
    --mod shell
    ```
    """
    def mapping_hla(self,sample):
        step = 'mapping_hla'
        cmd_line = self.get_cmd_line(step,sample)
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd = (
            f'{cmd_line} '
            f'--fq None '
        )
        self.process_cmd(cmd, step, sample, m=2, x=self.args.thread)


def main():
    # TODO
    multi = Multi_hla(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
