import os
import json
import io
import pandas as pd
from collections import abc
from jinja2 import Environment, PackageLoader, select_autoescape, FileSystemLoader
from celescope.tools.utils import add_log

class Reporter:
    def __init__(
            self, assay, step, sample, outdir, *,
            json_file=None, stat_file=None):

        self.assay = assay
        self.step = step
        self.sample = sample
        self.outdir = outdir
        if json_file:
            self.json_file = json_file
        else:
            self.json_file = f'{outdir}/../.metrics.json'

        if stat_file:
            self.stat_file = stat_file
        else:
            self.stat_file = f'{outdir}/stat.txt'

        if not os.path.exists(self.json_file):
            self.data = {}
        else:
            with open(self.json_file) as f:
                self.data = json.load(f)
        self.env = Environment(
            loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
            autoescape=select_autoescape(['html', 'xml'])
        )

    @add_log
    def dump_json(self):
        with open(self.json_file, 'w') as f:
            json.dump(self.data, f, indent=4)

    @add_log
    def render_html(self):
        template = self.env.get_template(f'html/{self.assay}/base.html')
        report_html = f"{self.outdir}/{self.sample}_report.html"
        with io.open(report_html, 'w', encoding='utf8') as f:
            html = template.render(self.data)
            f.write(html)

    def stat_to_json(self):
        self.data[self.step + '_summary'] = self.stat_to_dict()

    def stat_to_dict(self):
        df = pd.read_table(self.stat_file, header=None, sep=':', dtype=str)
        dic = dict(zip(df.iloc[:, 0], df.iloc[:, 1].str.strip()))
        metrics = dict()
        for key, value in dic.items():
            bool_fraction = False
            if '%' in value:
                bool_fraction = True
            chars = [',', '%', ')']
            for c in chars:
                value = value.replace(c, '')
            fraction = None
            try:
                number, fraction = value.split('(')
            except ValueError:
                number = value
            if not number.isnumeric():
                try:
                    number = float(number)
                    if bool_fraction:
                        number = round(number / 100, 4)
                except ValueError:
                    pass
            else:
                number = int(number)
            metrics[key] = number
            if fraction:
                fraction = round(float(fraction) / 100, 4)
                metrics[key + ' Fraction'] = fraction
        return metrics

    def add_data_item(self, **kwargs):
        for key, value in kwargs.items():
            # if value is a dict, and some value in this dict is float, format these value
            if isinstance(value, abc.Mapping):
                for value_key, value_value in value.items():
                    if isinstance(value_value, float):
                        value[value_key] = round(value_value, 2)

            self.data[key] = value




