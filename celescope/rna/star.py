import subprocess
import re
import pandas as pd

import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH
from celescope.tools.star_mixin import StarMixin, get_opts_star_mixin
from celescope.tools.step import Step 


class Star_rna(Step, StarMixin):
    """
    Features
    - Align R2 reads to the reference genome with STAR.
    - Collect Metrics with Picard.

    Output
    - `{sample}_Aligned.sortedByCoord.out.bam` BAM file contains Uniquely Mapped Reads.

    - `{sample}_SJ.out.tab` SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format.

    - `{sample}_Log.out` Main log with a lot of detailed information about the run. 
    This is most useful for troubleshooting and debugging.

    - `{sample}_Log.progress.out` Report job progress statistics, such as the number of processed reads, 
    % of mapped reads etc. It is updated in 1 minute intervals.

    - `{sample}_Log.Log.final.out` Summary mapping statistics after mapping job is complete, 
    very useful for quality control. The statistics are calculated for each read (single- or paired-end) and 
    then summed or averaged over all reads. Note that STAR counts a paired-end read as one read, 
    (unlike the samtools agstat/idxstats, which count each mate separately). 
    Most of the information is collected about the UNIQUE mappers 
    (unlike samtools agstat/idxstats which does not separate unique or multi-mappers). 
    Each splicing is counted in the numbers of splices, which would correspond to 
    summing the counts in SJ.out.tab. The mismatch/indel error rates are calculated on a per base basis, 
    i.e. as total number of mismatches/indels in all unique mappers divided by the total number of mapped bases.

    - `{sample}_region.log` Picard CollectRnaSeqMetrics results.
    """

    def __init__(self, args,display_title=None):
        Step.__init__(self, args, display_title=display_title)
        StarMixin.__init__(self, args)
        # parse
        self.refflat = f"{self.genomeDir}/{self.genome['refflat']}"

        self.ribo_log = f'{self.out_prefix}_ribo_log.txt'
        self.ribo_run_log = f'{self.out_prefix}_ribo_run.log'
        self.picard_region_log = f'{self.out_prefix}_region.log'
        self.plot = None
        self.stats = pd.Series()

    def add_other_metrics(self):
        """
        add picard region bases
        add region plot
        if debug, add ribosomal RNA reads percent
        """
        with open(self.STAR_map_log, 'r') as map_log:
            # number amd percent
            unique_reads_list = []
            multi_reads_list = []
            total_reads = 0
            for line in map_log:
                if line.strip() == '':
                    continue
                if re.search(r'Uniquely mapped reads', line):
                    unique_reads_list.append(line.strip().split()[-1])
                if re.search(r'of reads mapped to too many loci', line):
                    multi_reads_list.append(line.strip().split()[-1])
                if re.search(r'Number of input reads', line):
                    total_reads = int(line.strip().split()[-1])

        unique_reads = int(unique_reads_list[0])
        multi_reads = int(multi_reads_list[0])

        self.add_metric(
            name='Genome',
            value=self.genome['genome_name']
        )
        self.add_metric(
            name=f'Uniquely Mapped {self.stat_prefix}',
            value=unique_reads,
            total=total_reads,
            help_info='reads that mapped uniquely to the genome'
        )
        self.add_metric(
            name=f'Multi-Mapped {self.stat_prefix}',
            value=multi_reads,
            total=total_reads,
            help_info='reads that mapped to multiple locations in the genome'
        )

        with open(self.picard_region_log, 'r') as picard_log:
            region_dict = {}
            for line in picard_log:
                if not line:
                    break
                if line.startswith('## METRICS CLASS'):
                    header = picard_log.readline().strip().split('\t')
                    data = picard_log.readline().strip().split('\t')
                    region_dict = dict(zip(header, data))
                    break

        total = float(region_dict['PF_ALIGNED_BASES'])
        exonic_regions = int(region_dict['UTR_BASES']) + \
            int(region_dict['CODING_BASES'])
        intronic_regions = int(region_dict['INTRONIC_BASES'])
        intergenic_regions = int(region_dict['INTERGENIC_BASES'])

        self.add_metric(
            name='Base Pairs Mapped to Exonic Regions',
            value=exonic_regions,
            total=total,
            help_info='bases in primary alignments that align to a coding base or a UTR base for some gene'
        )
        self.add_metric(
            name='Base Pairs Mapped to Intronic Regions',
            value=intronic_regions,
            total=total,
            help_info='bases in primary alignments that align to an intronic base for some gene, and not a coding or UTR base'
        )
        self.add_metric(
            name='Base Pairs Mapped to Intergenic Regions',
            value=intergenic_regions,
            total=total,
            help_info='bases in primary alignments that do not align to any gene'
        )

        # ribo
        if self.debug:
            with open(self.ribo_log, 'r') as ribo_log:
                for line in ribo_log:
                    if line.find('#Matched') != -1:
                        items = line.split()
                        Reads_Mapped_to_rRNA = int(items[1])
                    if line.find('#Total') != -1:
                        items = line.split()
                        Reads_Total = int(items[1])
                self.add_metric(
                    name=f'{self.stat_prefix} Mapped to rRNA',
                    value=Reads_Mapped_to_rRNA,
                    total=Reads_Total,
                    help_info='Number of reads or umis that mapped to rRNA'
                )

        region_plot = {'region_labels': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
                       'region_values': [exonic_regions, intronic_regions, intergenic_regions]}
        self.add_data(region_plot=region_plot)

    @utils.add_log
    def ribo(self):
        # TODO remove bbduk.sh and use picard ribo bases
        human_ribo_fa = f'{ROOT_PATH}/data/rRNA/human_ribo.fasta'
        cmd = (
            f'bbduk.sh '
            f'in1={self.fq} '
            f'ref={human_ribo_fa} '
            f'stats={self.ribo_log} '
            f'overwrite=t '
            f'> {self.ribo_run_log} 2>&1 '
        )
        Star_rna.ribo.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def picard(self):
        cmd = [
            'picard',
            '-Xmx20G',
            '-XX:ParallelGCThreads=4',
            'CollectRnaSeqMetrics',
            'I=%s' % (self.STAR_bam),
            'O=%s' % (self.picard_region_log),
            'REF_FLAT=%s' % (self.refflat),
            'STRAND=NONE',
            'VALIDATION_STRINGENCY=SILENT']
        cmd_str = ' '.join(cmd)
        Star_rna.picard.logger.info(cmd_str)
        subprocess.check_call(cmd)

    @utils.add_log
    def run(self):
        self.run_star()
        self.picard()
        if self.debug:
            self.ribo()
        self.add_other_metrics()


def star(args):
    with Star_rna(args,display_title="Mapping") as runner:
        runner.run()

def get_opts_star(parser, sub_program):
    get_opts_star_mixin(parser, sub_program)
