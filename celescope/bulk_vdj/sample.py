from celescope.tools import utils
from celescope.__init__ import __VERSION__
from celescope.tools.step import Step, s_common

from xopen import xopen
import pysam 


class Sample(Step):
    """

    ## Features
    - Generate sample info.
    - Add read_ID in fastq file. The format of the read name is `{readId}_{index}`.

    ## Output
    - `01.barcode/{sample}_2.fq(.gz)`

    """
    def __init__(self, args):
        Step.__init__(self, args)
        self.assay_description = "Bulk-vdj"
        self.version = __VERSION__
        self.chemistry = "Customized"

        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")
        assert len(self.fq1_list) == len(self.fq2_list)
        self.fq_number = len(self.fq2_list)
        self.total_num = 0
        self.output_R1 = args.output_R1

        if args.gzip:
            self.suffix = ".gz"
        else:
            self.suffix = ""
        self.out_fq2 = f'{self.out_prefix}_2.fq{self.suffix}'
        #self.out_fq2 = f'{self.out_prefix}_clean_2.fa{self.suffix}'
        self.out_fq1 = f'{self.out_prefix}_1.fq{self.suffix}'

        self.open_files()

    def open_files(self):
        if self.output_R1:
            self.fh_fq1 = xopen(self.out_fq1, 'w')
        self.fh_fq2 = xopen(self.out_fq2, 'w')

    def close_files(self):
        if self.output_R1:
            self.fh_fq1.close()
        self.fh_fq2.close()


    @utils.add_log
    def add_readID(self):
        for i in range(self.fq_number):
            with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1, \
                    pysam.FastxFile(self.fq2_list[i], persist=False) as fq2:
                for entry1, entry2 in zip(fq1, fq2):
                    _, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                    _, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    self.total_num += 1

                    self.fh_fq2.write(f'@readId_{self.total_num}\n{seq2}\n+\n{qual2}\n')
                    #self.fh_fq2.write(f'>readId_{self.total_num}\n{seq1}{seq2}\n')
                    if self.output_R1:
                        self.fh_fq1.write(f'@readId_{self.total_num}\n{seq1}\n+\n{qual1}\n')       
                
                self.run.logger.info(self.fq1_list[i] + ' finished.')

        self.close_files()

    @utils.add_log
    def run(self):

        self.add_metric(
            name='Sample ID',
            value=self.sample,
        )
        self.add_metric(
            name='Assay',
            value=self.assay_description,
        )
        self.add_metric(
            name="Chemistry",
            value=self.chemistry,
        )
        self.add_metric(
            name="Software Version",
            value=self.version,
        )

        self.add_readID()

        self.add_metric(
            name = "Raw reads",
            value =self.total_num,
        )


@utils.add_log
def sample(args):
    with Sample(args) as runner:
        runner.run()


def get_opts_sample(parser, sub_program):
    parser.add_argument('--gzip', help="Output gzipped fastq files.", action='store_true')
    parser.add_argument('--output_R1', help="Output valid R1 reads.", action='store_true')
    if sub_program:
        parser.add_argument('--fq1', help='R1 fastq file. Multiple files are separated by comma.', required=True)
        parser.add_argument('--fq2', help='R2 fastq file. Multiple files are separated by comma.', required=True)
        parser = s_common(parser)
    return parser
