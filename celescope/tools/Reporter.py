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
            stat_file=None):

        self.assay = assay
        self.step = step
        self.sample = sample
        self.outdir = outdir
        self.path_dict = {
            "metric": f'{outdir}/../.metrics.json',
            "data": f'{outdir}/../.data.json'
        }

        if stat_file:
            self.stat_file = stat_file
        else:
            self.stat_file = f'{outdir}/stat.txt'

        # check dir
        if not os.path.exists(self.outdir):
            os.system('mkdir -p %s' % self.outdir)

        self.content_dict = {}
        for slot, path in self.path_dict.items():
            if not os.path.exists(path):
                self.content_dict[slot] = {}
            else:
                with open(path) as f:
                    self.content_dict[slot] = json.load(f)

        self.env = Environment(
            loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
            autoescape=select_autoescape(['html', 'xml'])
        )

    def dump_content(self, slot):
        '''dump content to json file
        '''
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
        self.content_dict['data'][self.step + '_summary'] = df.values.tolist()

    def stat_to_metric(self):
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
        self.content_dict['metric'][self.step + '_summary'] = metrics

    def add_content_item(self, slot, **kwargs):
        for key, value in kwargs.items():
            # if value is a dict, and some value in this dict is float, format these value
            if isinstance(value, abc.Mapping):
                for value_key, value_value in value.items():
                    if isinstance(value_value, float):
                        value[value_key] = round(value_value, 4)

            self.content_dict[slot][key] = value

    def add_data_item(self, **kwargs):
        self.add_content_item("data", **kwargs)

    @staticmethod
    def get_table(title, id, df_table):
        """
        return html code
        """
        table_dict = {}
        table_dict['title'] = title
        table_dict['table'] = df_table.to_html(
            escape=False,
            index=False,
            table_id=id,
            justify="center")
        table_dict['id'] = id
        return table_dict




