from celescope.tools.multi import Multi
from celescope.vdj.__init__ import __ASSAY__


class Multi_vdj(Multi):
    """

    ## Reference
    - Human
    ```
    mkdir -p /genome/vdj/human
    cd /genome/vdj/human
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TR{A,B}{V,J}.fasta
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IG{H,K,L}{V,J}.fasta
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta
    celescope vdj mkref human TR
    celescope vdj mkref human IG
    ```

    - Mouse
    ```
    mkdir -p /genome/vdj/mouse
    cd /genome/vdj/mouse
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TR{A,B}{V,J}.fasta
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/TRBD.fasta
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IG{H,K,L}{V,J}.fasta
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGHD.fasta
    celescope vdj mkref mouse TR
    celescope vdj mkref mouse IG
    ```
 
    ## Usage

    Please note `multi_vdj` is only used for CDR3 vdj data(not full length). For full length vdj data, please use `multi_flv_CR` or  `multi_flv_trust4`.
    
    ```
    multi_vdj \\
        --mapfile ./vdj.mapfile \\
        --ref_path /genome/vdj/human/human_TR \\
        --species human \\
        --type TCR \\
        --thread 8 \\
        --mod shell
    ``` 

    """

    def consensus(self, sample):
        step = "consensus"
        fq = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = f"{cmd_line} " f"--fq {fq} " f"--out_fasta "
        self.process_cmd(cmd, step, sample, m=5, x=1)
        outfile = f"{self.outdir_dic[sample][step]}/{sample}_consensus.fasta"
        return outfile

    def mapping_vdj(self, sample):
        step = "mapping_vdj"
        cmd_line = self.get_cmd_line(step, sample)
        fasta = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fasta'
        cmd = f"{cmd_line} " f"--fasta {fasta} "
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)

    def count_vdj(self, sample):
        # count_vdj
        step = "count_vdj"
        cmd_line = self.get_cmd_line(step, sample)
        UMI_count_filter_file = (
            f'{self.outdir_dic[sample]["mapping_vdj"]}/{sample}_UMI_count_filtered.tsv'
        )
        cmd = (
            f"{cmd_line} "
            f"--match_dir {self.col4_dict[sample]} "
            f"--UMI_count_filter_file {UMI_count_filter_file} "
        )
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)


def main():
    multi = Multi_vdj(__ASSAY__, min_col=4)
    multi.run()


if __name__ == "__main__":
    main()
