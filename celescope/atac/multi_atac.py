from celescope.tools.multi import Multi, TOOLS_DIR ,get_read
from celescope.atac.__init__ import __ASSAY__
from celescope.atac.atac import __SUB_STEPS__
from collections import defaultdict
from celescope.tools import utils
from itertools import combinations


def get_fq(library_id, library_path):
    fq1_list = get_read(library_id, library_path, read='1')
    fq2_list = get_read(library_id, library_path, read='2')
    fq3_list = get_read(library_id, library_path, read='3')
    for fq_list in combinations([fq1_list, fq2_list, fq3_list], 2):
        if len(fq_list[0]) != len(fq_list[1]):
            raise Exception("Read fastq number do not match!")

    fq1 = ",".join(fq1_list)
    fq2 = ",".join(fq2_list)
    fq3 = ",".join(fq3_list)
    return fq1, fq2, fq3


class Multi_atac(Multi):
    """
    ## Download and unpack cellranger soft and reference file.
    ```
    wget -O cellranger-atac-2.1.0.tar.gz "https://cf.10xgenomics.com/releases/cell-atac/cellranger-atac-2.1.0.tar.gz?Expires=1659105655&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1hdGFjL2NlbGxyYW5nZXItYXRhYy0yLjEuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTkxMDU2NTV9fX1dfQ__&Signature=ZzSn3fB9aI4GR9S4uaUOSaH3aKV5aJG2JgPTIQIV5Uka3GYhjc8QSTU~Eb4osKlLo8pghC4ze0PqpwOwUW6UGbhaX~eSoj-Vvei~kSxX3f0psyJ6Kbi4yscWe48RFIRaZkL91pFdApjEskhvDVh2a5eKTpBkt3dolSLncmyoc~7JT-D3dDymYvkZiM2KcfPHBi2lkBp5kxsYVmEzwH3btE0EczIjPCJwZbXHn7CU11XxyrTkp0gUU4Yp2QzglzggJ1kPveSlaUxeZgv~YVs9d-VQ36yVhgu~Yum~tcDZ~1hyr1aR9uigllni8g90uTdKfEaiC6idJ3kh2k9mlWh6iQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
    Reference: human, mouse, human-and-mouse
    wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
    wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
    wget https://cf.10xgenomics.com/supp/cell-atac/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0.tar.gz
    tar -xzvf cellranger-atac-2.1.0.tar.gz
    tar -xzvf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
    tar -xzvf refdata-cellranger-arc-mm10-2020-A-2.0.0.tar.gz
    tar -xzvf refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0.tar.gz
    ```
    ## Usage
    
    ```
    conda activate celescope
    multi_atac \\
        --mapfile ./test.mapfile \\
        --customized \\
        --pattern \\
        --whitelist \\
        --linker \\
        --ref_path "/SGRNJ06/randd/USER/cjj/cr-atac/refdata-cellranger-arc-GRCh38-2020-A-2.0.0" \\
        --soft_path "/SGRNJ06/randd/USER/cjj/cr-atac/cellranger-atac-2.1.0/" \\
        --mod shell
    ```
    """

    @staticmethod
    @utils.add_log
    def parse_mapfile(mapfile, default_val):
        fq_dict = defaultdict(list)
        col4_dict = {}
        col5_dict = {}
        with open(mapfile) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                line_split = line.split()
                library_id, library_path, sample_name = line_split[:3]
                if len(line_split) >= 4:
                    col4 = line_split[3]
                else:
                    col4 = default_val
                fq1, fq2, fq3 = get_fq(library_id, library_path)

                if sample_name in fq_dict:
                    fq_dict[sample_name][0].append(fq1)
                    fq_dict[sample_name][1].append(fq2)
                    fq_dict[sample_name][2].append(fq3)
                else:
                    fq_dict[sample_name] = [[fq1], [fq2], [fq3]]
                    col4_dict[sample_name] = col4
                if len(line_split) == 5:
                    col5_dict[sample_name] = line_split[4]

        for sample_name in fq_dict:
            fq_dict[sample_name][0] = ",".join(fq_dict[sample_name][0])
            fq_dict[sample_name][1] = ",".join(fq_dict[sample_name][1])
            fq_dict[sample_name][2] = ",".join(fq_dict[sample_name][2])

        if not fq_dict:
            raise Exception('empty mapfile!')
        return fq_dict, col4_dict, col5_dict

    def barcode(self, sample):
        step = "barcode"
        arr = self.fq_dict[sample]
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {arr[0]} --fq2 {arr[1]} --fq3 {arr[2]}'
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)
        
    def convert(self, sample):
        step = 'convert'
        fq1 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_1.fq'
        fq2 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_2.fq'
        fq3 = f'{self.outdir_dic[sample]["barcode"]}/{sample}_3.fq'
        cmd_line = self.get_cmd_line(step, sample)
        cmd = (
            f'{cmd_line} '
            f'--fq1 {fq1} --fq2 {fq2} --fq3 {fq3} '
        )
        self.process_cmd(cmd, step, sample, m=5, x=1)
    
    def atac(self, sample):
        step = 'atac'
        cmd_line = self.get_cmd_line(step, sample)
        fqs_dir = f'{self.outdir_dic[sample]["convert"]}'
        cmd = (
            f'{cmd_line} '
            f'--fqs_dir {fqs_dir} '
        )
        self.process_cmd(cmd, step, sample, m=self.args.mem, x=self.args.thread)

    def merge_report(self):
        step = "merge_report"
        _index = self.__STEPS__.index('atac') + 1
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
    multi = Multi_atac(__ASSAY__)
    multi.run()
    

if __name__ == '__main__':
    main()