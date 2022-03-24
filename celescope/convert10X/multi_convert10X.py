from celescope.tools.multi import Multi
from celescope.convert10X.__init__ import __ASSAY__


class Multi_convert10X(Multi):
    """

    ## Installation

    1. Clone repo
    ```
    git clone -b convert_10X https://github.com/singleron-RD/CeleScope.git
    ```

    2. Create conda environment and install conda packages
    ```
    cd CeleScope
    conda create -n convert_10X -y --file conda_pkgs.txt
    ```

    Alternatively, you can use [mamba](https://github.com/mamba-org/mamba) to improve speed.
    ```
    conda install mamba
    mamba create -n convert_10X -y --file conda_pkgs.txt
    ```

    3. Install celescope

    Make sure you have activated the conda environment before running `pip install Celescope`. 
    ```
    conda activate convert_10X
    pip install .
    ```

    ## Usage
    
    ```
    conda activate convert_10X
    multi_convert10X \\
        --mapfile  rna.mapfile \\
        --chemistry flv_rna \\
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


def main():
    multi = Multi_convert10X(__ASSAY__)
    multi.run()
    

if __name__ == '__main__':
    main()
