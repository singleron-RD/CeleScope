import io
import json
import os

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

env = Environment(
    loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
    autoescape=select_autoescape(['html', 'xml']),
)


class reporter:

    def __init__(
        self, assay, name, outdir, sample,
        stat_file=None, plot=None, parameters=None,
        table_file=None, table_header=None, df=None, html_flag=True
    ):
        self.assay = assay
        self.name = name
        self.stat_file = stat_file
        self.outdir = outdir
        self.sample = sample
        self.plot = plot
        self.parameters = parameters
        self.table_file = table_file
        self.table_header = table_header
        self.df = df
        self.html_flag = html_flag

    def get_report(self):

        json_file = self.outdir + '/.data.json'
        if not os.path.exists(json_file):
            data = {}
        else:
            fh = open(json_file)
            data = json.load(fh)
            fh.close()

        if self.stat_file:
            df = pd.read_table(self.stat_file, header=None, sep=':', dtype=str)
            data[self.name + '_summary'] = df.values.tolist()

        if self.plot:
            data[self.name + '_plot'] = self.plot

        if self.parameters:
            data[self.name + '_parameters'] = self.parameters

        if self.table_file:
            df = pd.read_csv(self.table_file, sep="\t")
            df = df.fillna(value="")
            data[self.name + '_table'] = df.values.tolist()

        if isinstance(self.df, pd.DataFrame):
            df = self.df.fillna(value="")
            data[self.name + '_table'] = df.values.tolist()

        if self.table_header:
            data[self.name + '_table_header'] = self.table_header

        if self.html_flag:
            template = env.get_template(f'html/{self.assay}/base.html')
            report_html = "{outdir}/{sample}_report.html".format(
                outdir=self.outdir, sample=self.sample)
            with io.open(report_html, 'w', encoding='utf8') as fh:
                html = template.render(data)
                # fh.write(html.encode('utf-8'))
                fh.write(html)

        with open(json_file, 'w') as fh:
            json.dump(data, fh, indent=4)
