import abc
import sys
import io
import json
import numbers
import os
import subprocess

from jinja2 import Environment, FileSystemLoader, select_autoescape

import celescope.tools.utils as utils
from celescope.__init__ import HELP_DICT


def cap_str_except_preposition(my_string):
    prepositions = {"and", "or", "the", "a", "of", "in", "per", "after", 'with'}
    lowercase_words = my_string.split(" ")

    final_words = [word if word in prepositions else word[0].upper() + word[1:] for word in lowercase_words]
    final_words = " ".join(final_words)
    return final_words


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

    def __init__(self, args, display_title=None):
        '''
        display_title controls the section title in HTML report
        '''

        self.args = args
        self.outdir = args.outdir
        self.sample = args.sample
        self.assay = args.assay
        self.thread = int(args.thread)
        self.debug = args.debug
        self.out_prefix = f'{self.outdir}/{self.sample}'

        # important! make outdir before path_dict because path_dict use relative path.
        utils.check_mkdir(self.outdir)

        # set
        class_name = self.__class__.__name__
        if not display_title:
            self._display_title = class_name
        else:
            self._display_title = display_title
        self._step_name = class_name[0].lower() + class_name[1:]
        self.__slots = ['data', 'metrics']
        self._step_summary_name = f'{self._step_name}_summary'

        self.__metric_list = []
        self.__help_content = []
        self._path_dict = {}
        for slot in self.__slots:
            self._path_dict[slot] = f'{self.outdir}/../.{slot}.json'

        self.__content_dict = {}
        for slot, path in self._path_dict.items():
            if not os.path.exists(path):
                self.__content_dict[slot] = {}
            else:
                with open(path) as f:
                    self.__content_dict[slot] = json.load(f)
            # clear step_summary
            self.__content_dict[slot][self._step_summary_name] = {}

        # jinja env
        self.env = Environment(
            loader=FileSystemLoader(os.path.dirname(__file__) + '/../templates/'),
            autoescape=select_autoescape(['html', 'xml'])
        )

        # out file
        self.__stat_file = f'{self.outdir}/stat.txt'

    def add_metric(self, name, value, total=None, help_info=None, display=None):
        '''
        add metric to metric_list
        
        Args
            display: controls how to display the metric in HTML report.
        '''

        name = cap_str_except_preposition(name)
        if help_info:
            help_info = help_info[0].upper() + help_info[1:]
            if help_info[-1] != '.':
                help_info += '.'
        if not display:
            if isinstance(value, numbers.Number):
                display = str(format(value, ','))
            else:
                display = value
        fraction = None
        if total:
            fraction = round(value / total * 100, 2)
            display += f'({fraction}%)'
        self.__metric_list.append(
            {
                "name": name,
                "value": value,
                "total": total,
                "fraction": fraction,
                "display": display,
                "help_info": help_info,
            }
        )

    def _write_stat(self):
        with open(self.__stat_file, 'w') as writer:
            for metric in self.__metric_list:
                name = metric['name']
                display = metric['display']

                line = f'{name}: {display}'
                writer.write(line + '\n')

    def _dump_content(self):
        '''dump content to json file
        '''
        for slot, path in self._path_dict.items():
            if self.__content_dict[slot]:
                with open(path, 'w') as f:
                    json.dump(self.__content_dict[slot], f, indent=4)

    @utils.add_log
    def _render_html(self):
        template = self.env.get_template(f'html/{self.assay}/base.html')
        report_html = f"{self.outdir}/../{self.sample}_report.html"
        with io.open(report_html, 'w', encoding='utf8') as f:
            html = template.render(self.__content_dict['data'])
            f.write(html)

    def _add_content_data(self):
        step_summary = {}
        step_summary['display_title'] = self._display_title
        step_summary['metric_list'] = self.__metric_list
        step_summary['help_content'] = self.__help_content
        self.__content_dict['data'][self._step_summary_name].update(step_summary)

    def _add_content_metric(self):
        metric_dict = dict()
        for metric in self.__metric_list:
            name = metric['name']
            value = metric['value']
            fraction = metric['fraction']
            metric_dict[name] = value
            if fraction:
                metric_dict[f'{name} Fraction'] = fraction

        self.__content_dict['metrics'][self._step_summary_name].update(metric_dict)

    def add_data(self, **kwargs):
        """
        add data(other than metrics) to self.content_dict['data']
        for example: add plots and tables
        """
        for key, value in kwargs.items():
            self.__content_dict['data'][self._step_summary_name][key] = value

    def add_help_content(self, name, content):
        """
        add help info before metrics' help_info
        """
        self.__help_content.append(
            {
                'name': name,
                'content': content
            }
        )

    def get_slot_key(self, slot, step_name, key):
        '''read slot from json file
        '''
        return self.__content_dict[slot][step_name + '_summary'][key]

    def get_table_dict(self, title, table_id, df_table):
        """
        table_dict {title: '', table_id: '', df_table: pd.DataFrame}
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
    def _clean_up(self):
        self._add_content_data()
        self._add_content_metric()
        self._write_stat()
        self._dump_content()
        self._render_html()

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
        self._clean_up()
