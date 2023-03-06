from celescope.flv_trust4.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_flv_trust4(Multi):
    """

    ## Usage
    
    ```
    multi_flv_trust4 \\
        --mapfile ./test.mapfile \\
        --ref hg38 \\
        --thread 8 \\
        --seqtype BCR \\
        --mod shell
    ```
    """
    def barcode(self, sample):
        step = "barcode"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)

    def mapping(self, sample):
        step = 'mapping'
        cmd_line = self.get_cmd_line(step, sample)
        match_fq1 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_1.fq'
        match_fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--match_fq1 {match_fq1} '
            f'--match_fq2 {match_fq2} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=self.args.thread)

    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        candidate_fq = f'{self.outdir_dic[sample]["mapping"]}/{sample}_bcrtcr.fq'
        cmd = (
            f'{cmd_line} '
            f'--candidate_fq {candidate_fq} '
        )
        self.process_cmd(cmd, step, sample, m=30, x=self.args.thread)
    
    def summarize(self, sample):
        step = 'summarize'
        cmd_line = self.get_cmd_line(step, sample)
        assemble_out = f'{self.outdir_dic[sample]["assemble"]}/assemble'
        fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        cmd = (
            f'{cmd_line} '
            f'--assemble_out {assemble_out} '
            f'--fq2 {fq2} '
            f'--match_dir {self.col4_dict[sample]} '
        )
        self.process_cmd(cmd, step, sample, m=15, x=1)

    def annotation(self,sample):
        step = 'annotation'
        cmd_line = self.get_cmd_line(step,sample)
        summarize_out = f'{self.outdir_dic[sample]["summarize"]}'
        match_dir = f'{self.col4_dict[sample]}'
        cmd = (
            f'{cmd_line} '
            f'--summarize_out {summarize_out} '
            f'--match_dir {match_dir} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)


def main():
    multi = Multi_flv_trust4(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()