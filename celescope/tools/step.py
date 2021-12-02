import abc
import sys
import io
import json
import numbers
import os
import subprocess
from collections import namedtuple

from jinja2 import Environment, FileSystemLoader, select_autoescape

import celescope.tools.utils as utils
from celescope.__init__ import HELP_DICT


def s_common(parser):
    """subparser common arguments
    """
    parser.add_argument('--outdir', help='Output diretory.', required=True)
    parser.add_argument('--assay', help='Assay name.', required=True)
    parser.add_argument('--sample', help='Sample name.', required=True)
    parser.add_argument('--thread', help=HELP_DICT['thread'], default=4)
    parser.add_argument('--debug', help=HELP_DICT['debug'], action='store_true')
    return parser

class ExtendEncoder(json.JSONEncoder):
    """
    convert numpy data types to python data types
    json does not recognize NumPy data types
    """
    def default(self, obj):
        if isinstance(obj, (np.int64,np.int32)):
            return int(obj)
        elif isinstance(obj,(np.float32,np.float64)):
            return float(obj)
        elif isinstance(obj,(np.ndarray,)):
            return obj.tolist()
        else:
            return json.JSONEncoder.default(self, obj)

class Step:
    """
    Step class
    """

    def __init__(self, args, display_title=None):
        class_name = self.__class__.__name__
        if not display_title:
            self.display_title = class_name
        else:
            self.display_title = display_title
        self.step_name = class_name.lower()
        self.args = args
        self.outdir = args.outdir
        self.sample = args.sample
        self.assay = args.assay
        self.thread = int(args.thread)
        self.debug = args.debug


        # important! make outdir before path_dict because path_dict use relative path.
        utils.check_mkdir(self.outdir)

        # set
        self.__slots__ = ['data', 'metric']
        self.out_prefix = f'{self.outdir}/{self.sample}'
        self.step_summary_name = f'{self.step_name}_summary'
        self.metric_list = []
        self.help_content = []
        self.path_dict = {}
        for slot in self.__slots__:
            self.path_dict[slot] = f'{self.outdir}/../.{slot}.json'

        self.content_dict = {}
        for slot, path in self.path_dict.items():
            if not os.path.exists(path):
                self.content_dict[slot] = {}
            else:
                with open(path) as f:
                    self.content_dict[slot] = json.load(f)
            # clear step_summary
            self.content_dict[slot][self.step_summary_name] = {}

        # jinja env
        self.env = Environment(
            loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
            autoescape=select_autoescape(['html', 'xml'])
        )

        # out file
        self.stat_file = f'{self.outdir}/stat.txt'

    def add_metric(self, name, value, total=None, help_info=None, display=None):
        '''add metric to metric_list
        '''
        if not display:
            if isinstance(value, numbers.Number):
                display = str(format(value, ','))
            else:
                display = value
        fraction = None
        if total:
            fraction = round(value / total * 100, 2)
            display += f'({fraction}%)'
        self.metric_list.append(
            {
                "name": name,
                "value": value,
                "total": total,
                "fraction": fraction,
                "display": display,
                "help_info": help_info,
            }
        )

    def write_stat(self):
        with open(self.stat_file, 'w') as writer:
            for metric in self.metric_list:
                name = metric['name']
                display = metric['display']

                line = f'{name}: {display}'
                writer.write(line + '\n')

    def dump_content(self):
        '''dump content to json file
        '''
        for slot, path in self.path_dict.items():            
            if self.content_dict[slot]:
                with open(path, 'w') as f:
                    json.dump(self.content_dict[slot], f, indent=4, cls=ExtendEncoder)

    @utils.add_log
    def render_html(self):
        template = self.env.get_template(f'html/{self.assay}/base.html') 
        report_html = f"{self.outdir}/../{self.sample}_report.html"
        with io.open(report_html, 'w', encoding='utf8') as f:
            html = template.render(self.content_dict['data'])
            f.write(html)

    def add_content_data(self):
        step_summary = {}
        step_summary['display_title'] = self.display_title
        step_summary['metric_list'] = self.metric_list
        step_summary['help_content'] = self.help_content
        self.content_dict['data'][self.step_summary_name].update(step_summary)

    def add_content_metric(self):
        metric_dict = dict()
        for metric in self.metric_list:
            name = metric['name']
            value = metric['value']
            fraction = metric['fraction']
            metric_dict[name] = value
            if fraction:
                metric_dict[f'{name} Fraction'] = fraction

        self.content_dict['metric'][self.step_summary_name].update(metric_dict)
    
    def add_data(self, **kwargs):
        """
        add data(other than metrics) to self.content_dict['data']
        for example: add plots and tables
        """
        for key, value in kwargs.items():
            self.content_dict['data'][self.step_summary_name][key] = value

    def add_help_content(self, name, content):
        """
        add help info before metrics' help_info
        """
        self.help_content.append(
            {
                'name': name,
                'content': content
            }
        )

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

    @utils.add_log
    def clean_up(self):
        self.add_content_data()
        self.add_content_metric()
        self.write_stat()
        self.dump_content()
        self.render_html() 

    @utils.add_log
    def debug_subprocess_call(self, cmd):
        '''
        debug subprocess call
        '''
        self.debug_subprocess_call.logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    @abc.abstractmethod
    def run(self):
        sys.exit('Please implement run() method.')

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.clean_up()

