from celescope.tools.multi import Multi
from celescope.vdj.__init__ import __ASSAY__


class Multi_vdj(Multi):
    """
    ## Download and unpack igblast soft.

    Soft: IgBLAST v1.9.0 or higher is required. \\
    mkdir -p ~/biosoft/igblast \\
    cd ~/biosoft/igblast \\
    wget -c https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.20.0-x64-linux.tar.gz \\
    tar -xzf ncbi-igblast-1.20.0-x64-linux.tar.gz \\

    ## Download IMGT Reference from IMGT(http://www.imgt.org/) and build index for igblast
    mkdir -p ~/biosoft/imgt_ref \\
    cd ~/biosoft/imgt_ref \\
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TR{A,B}{V,J}.fasta \\
    wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta \\

    Combine all V, all D and all J sequences, respectively, into separate files: \\
    cat TRAV.fasta TRBV.fasta > TRV.fasta \\
    cat TRAJ.fasta TRBJ.fasta > TRJ.fasta \\
    cat TRBD.fasta > TRD.fasta \\

    Build Index for igblast: \\
    perl ~biosoft/igblast/ncbi-igblast-1.20.0/bin/edit_imgt_file.pl TRV.fasta > TRV.fa \\
    ~biosoft/igblast/ncbi-igblast-1.20.0/bin/makeblastdb -parse_seqids -dbtype nucl -in TRV.fa \\

    perl ~biosoft/igblast/ncbi-igblast-1.20.0/bin/edit_imgt_file.pl TRD.fasta > TRD.fa \\
    ~biosoft/igblast/ncbi-igblast-1.20.0/bin/makeblastdb -parse_seqids -dbtype nucl -in TRD.fa \\

    perl ~biosoft/igblast/ncbi-igblast-1.20.0/bin/edit_imgt_file.pl TRBJ.fasta > TRJ.fa \\
    ~biosoft/igblast/ncbi-igblast-1.20.0/bin/makeblastdb -parse_seqids -dbtype nucl -in human_TRJ.fa \\

    ## Usage
    ```
    multi_vdj \\
        --mapfile ./vdj.mapfile \\
        --soft_path /SGRNJ06/randd/USER/cjj/igblast/ncbi-igblast-1.20.0 \\
        --ref_path /SGRNJ06/randd/USER/cjj/igblast/igblast_ref/hs_TR \\
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