from celescope.tools.multi import Multi
from celescope.bulk_vdj.__init__ import __ASSAY__


class Multi_bulk_vdj(Multi):
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
    ```
    multi_bulk_vdj \\
        --mapfile ./vdj.mapfile \\
        --ref_path /genome/vdj/human/human_TR \\
        --species human \\
        --type TCR \\
        --min_consensus_read 2 \\
        --allowNoLinker \\
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
        # fasta = f'{self.outdir_dic[sample]["consensus"]}/{sample}_consensus.fasta'
        fasta = (
            f'{self.outdir_dic[sample]["consensus"]}/{sample}_filtered_consensus.fasta'
        )
        metrics = f'{self.outdir_dic[sample]["consensus"]}/{sample}_metrics.tsv'
        cmd = f"{cmd_line} " f"--fasta {fasta} " f"--consensus_metrics_file {metrics}"
        self.process_cmd(cmd, step, sample, m=15, x=self.args.thread)


def main():
    multi = Multi_bulk_vdj(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
