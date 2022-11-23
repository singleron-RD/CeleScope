from celescope.tools.multi import Multi
from celescope.vdj.__init__ import __ASSAY__


class Multi_vdj(Multi):
    """

    ## Download IMGT Reference from IMGT(http://www.imgt.org/)
    USE Human TR IMGT As Example:
    ```
    mkdir -p ~/biosoft/imgt_ref \\
    cd ~/biosoft/imgt_ref \\
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TR{A,B}{V,J}.fasta \\
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta
    ```

    ## Build Index for IMGT_ref
    Make sure running this command in imgt_ref directory which contains all V,D,J of TRA/TRB or IGH/IGK/IGL reference downloaded from IMGT website.
    ```
    celescope vdj mkref human TR.

    ~/biosoft/imgt_ref/human_TR will be generated.
    ```

    ## Usage
    ```
    multi_vdj \\
        --mapfile ./vdj.mapfile \\
        --ref_path ~/biosoft/imgt_ref/human_TR \\
        --species human \\
        --type TCR \\
        --thread 8 \\
        --mod shell
    ``` 

    """

    def consensus(self, sample):
        step = 'consensus'
        fq = f'{self.outdir_dic[sample]["cutadapt"]}/{sample}_clean_2.fq{self.fq_suffix}'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq {fq} '
            f'--out_fasta '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)
        outfile = f'{self.outdir_dic[sample][step]}/{sample}_consensus.fasta'
        return outfile

    def mapping_vdj(self, sample):
        step = 'mapping_vdj'
        cmd_line = self.get_cmd_line(step, sample)
        fasta = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fasta'
        cmd = (
            f'{cmd_line} '
            f'--fasta {fasta} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)

    def count_vdj(self, sample):
        # count_vdj
        step = 'count_vdj'
        cmd_line = self.get_cmd_line(step, sample)
        UMI_count_filter_file = (
            f'{self.outdir_dic[sample]["mapping_vdj"]}/{sample}_UMI_count_filtered.tsv'
        )
        cmd = (
            f'{cmd_line} '
            f'--match_dir {self.col4_dict[sample]} '
            f'--UMI_count_filter_file {UMI_count_filter_file} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)


def main():
    multi = Multi_vdj(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()