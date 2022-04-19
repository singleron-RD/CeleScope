from celescope.fl_vdj_TRUST4_split.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_fl_vdj_TRUST4_split(Multi):
    """

    ## Installation

    1. Clone repo
    ```
    git clone https://github.com/singleron-RD/CeleScope.git
    ```

    2. Create conda environment and install conda packages
    ```
    cd CeleScope
    conda create -n celescope -y --file conda_pkgs.txt
    ```

    Alternatively, you can use [mamba](https://github.com/mamba-org/mamba) to improve speed.
    ```
    conda install mamba
    mamba create -n celescope -y --file conda_pkgs.txt
    ```

    3. Install celescope

    Make sure you have activated the conda environment before running `pip install Celescope`. 
    ```
    conda activate celescope
    pip install .
    ```

    ## Usage
    
    ```
    conda activate TRUST_dev
        multi_fl_vdj_TRUST4_split \\
        --mapfile ./test.mapfile \\
        --outdir ./ \\
        --chemistry flv \\
        --allowNoLinker \\
        --species GRCm38 \\
        --thread 8 \\
        --seqtype BCR \\
    ```
    """
    
    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'

        cmd = (
            f'{cmd_line} '
            f'--fq2 {fq2} '
            f'--match_dir {self.col4_dict[sample]}'
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)
    
    def summarize(self, sample):
        step = 'summarize'
        cmd_line = self.get_cmd_line(step, sample)
        reads_assignment = f'{self.outdir_dic[sample]["assemble"]}/assemble/{sample}_assign.out'
        assembled_fa = f'{self.outdir_dic[sample]["assemble"]}/assemble/{sample}_assembled_reads.fa'
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--reads_assignment {reads_assignment} '
            f'--assembled_fa {assembled_fa} '
            f'--fq2 {fq2} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)
        
    def mapping_annotation(self,sample):
        step = 'mapping_annotation'
        cmd_line = self.get_cmd_line(step,sample)
        match_dir = f'{self.col4_dict[sample]}'
        cmd=(
            f'{cmd_line} '
            f'--match_dir {match_dir} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_fl_vdj_TRUST4_split(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
