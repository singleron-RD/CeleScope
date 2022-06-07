from celescope.fl_vdj_TRUST4.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_fl_vdj_TRUST4(Multi):
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
        multi_fl_vdj_TRUST4 \\
        --mapfile ./test.mapfile \\
        --outdir ./ \\
        --chemistry flv \\
        --allowNoLinker \\
        --species GRCm38 \\
        --thread 10 \\
        --seqtype BCR \\
        --match_previous_assemble
    ```
    """
    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'

        cmd = (
            f'{cmd_line} '
            f'--cutadapted_fq {fq2} '
            f'--match_dir {self.col4_dict[sample]}'
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)
    
    def summarize(self, sample):
        step = 'summarize'
        cmd_line = self.get_cmd_line(step, sample)
        full_len_assembly = f'{self.outdir_dic[sample]["assemble"]}/{sample}_full_len.fa'
        fq2 = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq'
        assign_out = f'{self.outdir_dic[sample]["assemble"]}/{sample}_assign.out'
        filter_report = f'{self.outdir_dic[sample]["assemble"]}/{sample}_filter_report.tsv'
        barcode_filter_report = f'{self.outdir_dic[sample]["assemble"]}/{sample}_barcode_filter_report.tsv'
        assembled_fa = f'{self.outdir_dic[sample]["assemble"]}/{sample}_assembled_reads.fa'

        cmd = (
            f'{cmd_line} '
            f'--full_len_assembly {full_len_assembly} '
            f'--cutadapted_fq {fq2} '
            f'--assign_out {assign_out} '
            f'--filter_report {filter_report} '
            f'--barcode_filter_report {barcode_filter_report} '
            f'--assembled_fa {assembled_fa} '
        )
        self.process_cmd(cmd, step, sample, m=10, x=5)

    def annotation(self,sample):
        step = 'annotation'
        cmd_line = self.get_cmd_line(step,sample)
        match_dir = f'{self.col4_dict[sample]}'
        cmd=(
            f'{cmd_line} '
            f'--match_dir {match_dir} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_fl_vdj_TRUST4(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()