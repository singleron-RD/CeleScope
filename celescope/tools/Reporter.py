import os
import json
import io
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape, FileSystemLoader
from celescope.tools.utils import add_log

class Reporter:
    def __init__(self, assay, step, sample, outdir, json_file):
        self.assay = assay
        self.step = step
        self.sample = sample
        self.outdir = outdir
        self.json_file = json_file

        if not os.path.exists(json_file):
            self.data = {}
        else:
            with open(json_file) as f:
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

    def stat_to_json(self, stat_file):
        df = pd.read_table(stat_file, header=None, sep=':', dtype=str)
        self.data[self.step + '_summary'] = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

    def add_data_item(self, **kwargs):
        for key, value in kwargs.items():
            self.data[key] = value




