import abc
import collections
import io
import json
import numbers
import os
from collections import namedtuple

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from celescope.tools.utils import add_log
from celescope.__init__ import HELP_DICT


Metric = namedtuple("Metric", "name value total fraction")


def s_common(parser):
    """subparser common arguments
    """
    parser.add_argument('--outdir', help='Output diretory.', required=True)
    parser.add_argument('--assay', help='Assay name.', required=True)
    parser.add_argument('--sample', help='Sample name.', required=True)
    parser.add_argument('--thread', help=HELP_DICT['thread'], default=4)
    parser.add_argument('--debug', help=HELP_DICT['debug'], action='store_true')
    return parser


class Step:
    """
    Step class
    """

    def __init__(self, args, step_name):
        self.step_name = step_name
        self.args = args
        self.outdir = args.outdir
        self.sample = args.sample
        self.assay = args.assay
        self.thread = int(args.thread)
        self.debug = args.debug
        # set
        self.out_prefix = f'{self.outdir}/{self.sample}'

        # important! make outdir before path_dict because path_dict use relative path.
        if not os.path.exists(self.outdir):
            os.system('mkdir -p %s' % self.outdir)

        self.metric_list = []
        self.path_dict = {
            "metric": f'{self.outdir}/../.metrics.json',
            "data": f'{self.outdir}/../.data.json'
        }
        self.content_dict = {}
        for slot, path in self.path_dict.items():
            if not os.path.exists(path):
                self.content_dict[slot] = {}
            else:
                with open(path) as f:
                    self.content_dict[slot] = json.load(f)

        # jinja env
        self.env = Environment(
            loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
            autoescape=select_autoescape(['html', 'xml'])
        )

        # out file
        self.stat_file = f'{self.outdir}/stat.txt'

    def add_metric(self, name, value=None, total=None, fraction=None):
        '''add metric to metric_list
        '''
        self.metric_list.append(Metric(
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
            metric_list.append(Metric(
                name=metric.name,
                value=metric.value,
                total=metric.total,
                fraction=fraction,
            ))
        self.metric_list = metric_list

    def metric_list_to_stat(self):
        with open(self.stat_file, 'w') as stat_handle:
            for metric in self.metric_list:
                line = f'{metric.name}: '
                value = metric.value
                fraction = metric.fraction
                value_bool = value == 0 or value
                fraction_bool = fraction == 0 or fraction
                if fraction_bool:
                    fraction = round(fraction * 100, 2)
                if value_bool:
                    if isinstance(value, numbers.Number):
                        line += format(value, ',')
                        if fraction_bool:
                            line += f'({fraction}%)'
                    else:
                        line += value
                elif fraction_bool:
                    line += f'{fraction}%'
                stat_handle.write(line + '\n')

    def dump_content(self, slot):
        '''dump content to json file
        '''
        if self.content_dict[slot]:
            with open(self.path_dict[slot], 'w') as f:
                json.dump(self.content_dict[slot], f, indent=4)

    @add_log
    def render_html(self):
        template = self.env.get_template(f'html/{self.assay}/base.html')
        report_html = f"{self.outdir}/../{self.sample}_report.html"
        with io.open(report_html, 'w', encoding='utf8') as f:
            html = template.render(self.content_dict['data'])
            f.write(html)

    def stat_to_data(self):
        df = pd.read_table(self.stat_file, header=None, sep=':', dtype=str)
        self.content_dict['data'][self.step_name + '_summary'] = df.values.tolist()

    def stat_to_metric(self):
        '''
        can be:
        1. value
        2. value(fraction%)
        3. fraction%
        '''

        df = pd.read_table(self.stat_file, header=None, sep=':', dtype=str)
        dic = dict(zip(df.iloc[:, 0], df.iloc[:, 1].str.strip()))
        metrics = dict()
        for metric_name, string in dic.items():
            bool_fraction = False
            bool_value = False
            if '%' in string:
                bool_fraction = True
                if "(" in string:
                    bool_value = True
            chars = [',', '%', ')']
            for character in chars:
                string = string.replace(character, '')

            if bool_fraction:
                if bool_value:  # case 2
                    value, fraction = string.split('(')
                    fraction = round(float(fraction) / 100, 4)
                    metrics[metric_name] = int(value)
                    metrics[metric_name + ' Fraction'] = fraction
                else:  # case 3
                    fraction = round(float(string) / 100, 4)
                    metrics[metric_name] = fraction
            else:  # case 1
                value = string
                if '.' in string:
                    try:
                        value = round(float(string), 4)
                    except ValueError:
                        pass
                else:
                    try:
                        value = int(string)
                    except ValueError:
                        pass
                metrics[metric_name] = value

        self.content_dict['metric'][self.step_name + '_summary'] = metrics

    def add_content_item(self, slot, **kwargs):
        for key, value in kwargs.items():
            # if value is a dict, and some value in this dict is float, format these value
            if isinstance(value, collections.abc.Mapping):
                for value_key, value_value in value.items():
                    if isinstance(value_value, float):
                        value[value_key] = round(value_value, 4)

            self.content_dict[slot][key] = value

    def add_data_item(self, **kwargs):
        self.add_content_item("data", **kwargs)

    @staticmethod
    def get_table(title, table_id, df_table):
        """
        return html code
        """
        table_dict = {}
        table_dict['title'] = title
        table_dict['table'] = df_table.to_html(
            escape=False,
            index=False,
            table_id=table_id,
            justify="center")
        table_dict['id'] = table_id
        return table_dict

    @add_log
    def clean_up(self):
        if self.metric_list:
            self.get_fraction()
            self.metric_list_to_stat()
        if os.path.exists(self.stat_file):
            self.stat_to_metric()
            self.stat_to_data()
        self.dump_content(slot="data")
        self.dump_content(slot="metric")
        self.render_html()

    @abc.abstractmethod
    def run(self):
        return
