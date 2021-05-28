import subprocess

import pandas as pd

from celescope.tools.step import Step
from celescope.tools.star_mixin import StarMixin, get_opts_star_mixin
import celescope.tools.utils as utils
from celescope.__init__ import ROOT_PATH


class Star_rna(Step, StarMixin):
    """
    star rna
    """
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
        StarMixin.__init__(self, args)
        # parse
        self.refflat = f"{self.genomeDir}/{self.genome['refflat']}"

        self.ribo_log = f'{self.outdir}/{self.sample}_ribo_log.txt'
        self.ribo_run_log = f'{self.outdir}/{self.sample}_ribo_run.log'
        self.picard_region_log = f'{self.outdir}/{self.sample}_region.log'
        self.plot = None
        self.stats = pd.Series()

    def add_other_metrics(self):
        """
        add picard region bases
        add region plot
        if debug, add ribosomal RNA reads percent
        """

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
        )
        self.add_metric(
            name='Base Pairs Mapped to Intronic Regions',
            value=intronic_regions,
            total=total,
        )
        self.add_metric(
            name='Base Pairs Mapped to Intergenic Regions',
            value=intergenic_regions, 
            total=total,
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
                )

        region_plot = {'region_labels': ['Exonic Regions', 'Intronic Regions', 'Intergenic Regions'],
                'region_values': [exonic_regions, intronic_regions, intergenic_regions]}   
        self.add_content_item("data", STAR_plot=region_plot)


    @utils.add_log
    def ribo(self):
        human_ribo_fa = f'{ROOT_PATH}/data/rRNA/human_ribo.fasta'
        self.ribo_log = f'{self.outdir}/{self.sample}_ribo_log.txt'
        self.ribo_run_log = f'{self.outdir}/{self.sample}_ribo_run.log'
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
        self.clean_up()


def star(args):
    step_name = "star"
    runner = Star_rna(args, step_name)
    runner.run()


def get_opts_star(parser, sub_program):
    get_opts_star_mixin(parser, sub_program)