from celescope.tools.multi import Multi
from celescope.convert10X.__init__ import __ASSAY__


class Multi_convert10X(Multi):
    """

    ## Usage
    Running cellranger-count for sgr sc-RNA data:
    ```
    conda activate celescope
    multi_convert10X \\
        --mapfile  test.mapfile \\
        --thread 8 \\
        --ref_path "/soft/cellranger/refdata-gex-GRCh38-2020-A" \\
        --soft_path "/soft/cellranger/cellranger-6.1.2/cellranger" \\
        --mod shell 
    ``` 
    Converting sgr data to 10X format:
    ```
    conda activate celescope
    multi_convert10X \\
        --mapfile  test.mapfile \\
        --thread 8 \\
        --soft_path "/soft/cellranger/cellranger-6.1.2/cellranger" \\
        --steps_run sample,barcode,convert \\
        --mod shell 
    ``` 
    """

    def convert(self, sample):
        step = 'convert'
        cmd_line = self.get_cmd_line(step, sample)
        fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--fq2 {fq2} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def cellranger(self, sample):
        step = 'cellranger'
        cmd_line = self.get_cmd_line(step, sample)
        fqs_dir = f'{self.outdir_dic[sample]["convert"]}'
        cmd = (
            f'{cmd_line} '
            f'--fqs_dir {fqs_dir} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.mem, x=self.args.thread)

        
def main():
    multi = Multi_convert10X(__ASSAY__)
    multi.run()
    

if __name__ == '__main__':
    main()
