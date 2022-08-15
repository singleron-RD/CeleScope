from celescope.tools.multi import Multi, TOOLS_DIR
from celescope.flv_CR.__init__ import __ASSAY__
from celescope.flv_CR.assemble import __SUB_STEPS__


class Multi_flv_CR(Multi):
    """

    ## Download and unpack cellranger soft and reference file.
    ```
    wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-vdj/cellranger-6.1.2.tar.gz?Expires=1646072261&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC12ZGovY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDYwNzIyNjF9fX1dfQ__&Signature=Z-2m906CV5Rb1snIAga-QDSXYSZ8cNqCj1EECGP4uloU3qH~uCMH42MHf4TNnDL2zAsKA7cXsCsQYz0A9yJdNh7dfRT8ohpuAzASFx5Pj-bkqfw4p2tql55IIaPN0zqxyUuyZ9sfKl5qTQX82LoVolRpiBUL8dF9nr~bA2P1gJZ~xg1QssS7icR5MmTzvKKS5NYkezG8vWaTiEdXU0nuKI2ciZSX5GOMeIRW-YYR7mJwHmBbTVxe0o-uBuUtqor0Y98jdIv8Z~dwMjujRjrEShdCGNixTSonGzeS2~9CXqWquCJIOolqFFkcFHgXkD7ZWNfSXWbTxuF57rCsub98pA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

    Reference: human and mouse
    wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
    wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz

    tar -xzvf cellranger-6.1.2.tar.gz
    tar -xzvf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
    tar -xzvf refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0.tar.gz
    ```

    ## Usage
    
    ```
    conda activate celescope
    multi_flv_CR \\
        --mapfile ./test.mapfile \\
        --thread 8 \\
        --seqtype TCR \\
        --ref_path "/soft/cellranger/vdj_ref/6.0.0/hs/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" \\
        --soft_path "/soft/cellranger/cellranger-6.1.2/cellranger" \\
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

    def assemble(self, sample):
        step = 'assemble'
        cmd_line = self.get_cmd_line(step, sample)
        fqs_dir = f'{self.outdir_dic[sample]["convert"]}'
        cmd = (
            f'{cmd_line} '
            f'--fqs_dir {fqs_dir} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.mem, x=self.args.thread)

    def summarize(self,sample):
        step = 'summarize'
        cmd_line = self.get_cmd_line(step, sample)
        barcode_convert_json = f'{self.outdir_dic[sample]["convert"]}/barcode_convert.json'
        assemble_out = f'{self.outdir_dic[sample]["assemble"]}/{sample}/outs'
        cmd = (
            f'{cmd_line} '
            f'--barcode_convert_json {barcode_convert_json} '
            f'--assemble_out {assemble_out} '
        )
        self.process_cmd(cmd, step, sample, m=8, x=self.args.thread)

    def match(self,sample):
        step = 'match'
        cmd_line = self.get_cmd_line(step, sample)
        summarize_out = f'{self.outdir_dic[sample]["summarize"]}'
        match_dir = f'{self.col4_dict[sample]}'
        cmd = (
            f'{cmd_line} '
            f'--match_dir {match_dir} '
            f'--summarize_out {summarize_out} '
        )
        self.process_cmd(cmd, step, sample, m=8, x= self.args.thread)
    
    def mapping(self,sample):
        step = 'mapping'
        cmd_line = self.get_cmd_line(step,sample)
        match_dir = f'{self.col4_dict[sample]}'
        match_out = f'{self.outdir_dic[sample]["match"]}'
        cmd = (
            f'{cmd_line} '
            f'--match_dir {match_dir} '
            f'--match_out {match_out} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)
    
    def merge_report(self):
        step = "merge_report"
        _index = self.__STEPS__.index('assemble') + 1
        steps_str = ",".join(self.__STEPS__[:_index] + __SUB_STEPS__ + self.__STEPS__[_index:-1])
        samples = ','.join(self.fq_dict.keys())
        app = TOOLS_DIR + '/merge_table.py'
        cmd = (
            f'python {app} --samples {samples} '
            f'--steps {steps_str} --outdir {self.args.outdir}'
        )
        if self.args.rm_files:
            cmd += ' --rm_files'
        self.generate_cmd(cmd, step, sample="")
        for sample in self.fq_dict:
            self.sjm_order += f'order {step} after {self.last_step}_{sample}\n'

def main():
    multi = Multi_flv_CR(__ASSAY__)
    multi.run()
    

if __name__ == '__main__':
    main()

