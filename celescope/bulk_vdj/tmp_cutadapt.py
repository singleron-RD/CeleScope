from celescope.tools import cutadapt as super_ct
from celescope.tools import utils

import subprocess


class Cutadapt(super_ct.Cutadapt):
    """
    ## Features
    - Trim adapters in R2 reads with cutadapt. Default adapters includes:
        - polyT=A{18}, 18 A bases. 
        - p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA, Illumina p5 adapter.
    ## Output
    - `cutadapt.log` Cutadapt output log file.
    - `{sample}_clean_2.fa` R2 reads file without adapters in fasta format.
    
    """
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.out_fa2 = f'{self.outdir}/{self.sample}_clean_2.fa'

    @utils.add_log
    def run(self):
        adapter_args_str = " ".join(['-a ' + adapter for adapter in self.adapter_args])
        cmd = (
            'cutadapt '
            f'{adapter_args_str} '
            f'-n {len(self.adapter_args)} '
            f'-j {self.thread} '
            f'-m {self.args.minimum_length} '
            f'--nextseq-trim={self.args.nextseq_trim} '
            f'--overlap {self.args.overlap} '
            f'-l {self.args.insert} '
            f'-o {self.out_fa2} '
            f'--fasta '
            f'{self.args.fq} '
        )
        self.run.logger.info(cmd)
        # need encoding argument to return str
        results = subprocess.run(
            cmd, stdout=subprocess.PIPE,
            encoding='utf-8', check=True, shell=True
        )
        cutadapt_log = results.stdout
        with open(self.cutadapt_log_file, 'w') as f:
            f.write(cutadapt_log)
        self.add_cutadapt_metrics(cutadapt_log)


@utils.add_log
def cutadapt(args):
    with Cutadapt(args, display_title="Trimming") as runner:
        runner.run()


def get_opts_cutadapt(parser, sub_program):
    super_ct.get_opts_cutadapt(parser, sub_program)