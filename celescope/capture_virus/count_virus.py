from celescope.tools.capture.count_bam import Count_bam, get_opts_count_bam


def get_opts_count_virus(parser, sub_program):
    get_opts_count_bam(parser, sub_program)


class Count_virus(Count_bam):
    """
    ## Features

    - Taking the bam file as input, count the number of UMIs and reads mapped to the viral genome.

    ## Output
    - {sample_raw_read_count.json} : barcode - UMI - raw_reads_count
    """



def count_virus(args):

    with Count_virus(args) as runner:
        runner.run()
