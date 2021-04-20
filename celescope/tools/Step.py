import os
import numbers
from collections import namedtuple
from celescope.tools.utils import add_log
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
        self.report = Reporter(
            self.assay,
            self.step_name,
            self.sample,
            self.outdir,
        )

        # check dir
        if not os.path.exists(self.outdir):
            os.system('mkdir -p %s' % self.outdir)

    def add_metric(self, name, value=None, total=None, fraction=None):
        '''add metric to metric_list
        '''
        self.metric_list.append(self.Metric(
            name=name, value=value, total=total, fraction=fraction
        ))

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

    def metric_list_to_stat(self):
        f_stat = open(self.stat_file, 'w')
        for metric in self.metric_list:
            line = f'{metric.name}: '
            value = metric.value
            fraction = metric.fraction
            if fraction:
                fraction = round(fraction * 100, 2)
            if value:
                if isinstance(value, numbers.Number):
                    line += format(value, ',')
                    if fraction:
                        line += f'({fraction}%)'
                else:
                    line += value
            elif fraction:
                line += f'{fraction}%'
            f_stat.write(line + '\n')
        f_stat.close()

    @add_log
    def clean_up(self):
        if self.metric_list:
            self.get_fraction()
            self.metric_list_to_stat()
        self.report.stat_to_metric()
        self.report.stat_to_data()
        self.report.dump_content(slot="data")
        self.report.dump_content(slot="metric")
        self.report.render_html()


