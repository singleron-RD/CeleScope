import os
import numbers
from collections import namedtuple
from celescope.tools.report import reporter
from celescope.tools.Reporter import Reporter


def s_common(parser):
    """subparser common arguments
    """
    parser.add_argument('--outdir', help='output dir', required=True)
    parser.add_argument('--assay', help='assay', required=True)
    parser.add_argument('--sample', help='sample name', required=True)
    parser.add_argument('--thread', default=4)
    parser.add_argument('--debug', help='debug', action='store_true')
    return parser


class Step:
    '''
    Step class
    '''
    def __init__(self, args, step_name):
        self.step_name = step_name
        self.args = args
        self.outdir = args.outdir
        self.sample = args.sample
        self.assay = args.assay
        self.thread = args.thread
        self.debug = args.debug
        self.stat_file = f'{self.outdir}/stat.txt'
        self.metric_list = []
        self.Metric = namedtuple("Metric", "name value total fraction")

        # check dir
        if not os.path.exists(self.outdir):
            os.system('mkdir -p %s' % self.outdir)

    def add_metric(self, name, value=None, total=None, fraction=None):
        self.metric_list.append(self.Metric(
            name=name, value=value, total=total, fraction=fraction
        ))


    def _get_report(self):
        report = reporter(
            name=self.step_name,
            assay=self.assay,
            sample=self.sample,
            stat_file=self.stat_file,
            outdir=self.outdir + '/..'
        )
        report.get_report()

    def get_report(self):
        '''
        need stat.txt
        '''
        report = Reporter(
            self.assay,
            self.step_name,
            self.sample,
            self.outdir,
            json_file = f'{self.outdir}/../.metrics.json'
        )
        report.stat_to_json()
        report.dump_json()

    def get_fraction(self):
        '''
        metric_list: list of namedtuple(name, value, total, fraction)
        '''
        metric_list = []
        for metric in self.metric_list:
            fraction = metric.fraction
            if metric.total:
                fraction = metric.value / metric.total
            if fraction:
                fraction = round(fraction, 4)
            metric_list.append(self.Metric(
                name=metric.name,
                value=metric.value,
                total=metric.total,
                fraction=fraction,
            ))
        self.metric_list = metric_list

    def get_stat(self):
        f_stat = open(self.stat_file, 'w')
        for metric in self.metric_list:
            line = f'{metric.name}: '
            value = metric.value
            fraction = metric.fraction
            if fraction:
                fraction = fraction * 100            
            if value:
                if isinstance(value, numbers.Number):
                    line += format(value, ',')
                    if fraction:
                        line += f'({metric.fraction}%)'
                else:
                    line += value
            elif fraction:
                line += f'{fraction}%'
            f_stat.write(line + '\n')
        f_stat.close()

    def clean_up(self):
        self.get_fraction()
        self.get_stat()
        self._get_report()
        self.get_report()

    #def tuple_to_metrcis
