import abc
from pathlib import Path
import shutil
import sys
import io
import json
import numbers
import os
import subprocess

from jinja2 import Environment, FileSystemLoader, select_autoescape

from celescope.tools import utils
from celescope.__init__ import HELP_DICT, __version__


def cap_str_except_preposition(my_string):
    prepositions = {
        "and",
        "or",
        "the",
        "a",
        "of",
        "in",
        "per",
        "after",
        "with",
        "by",
        "at",
    }
    lowercase_words = my_string.split(" ")

    final_words = [
        word if word in prepositions else word[0].upper() + word[1:]
        for word in lowercase_words
    ]
    final_words = " ".join(final_words)
    return final_words


def safe_move(source, dest_dir):
    source = Path(source)
    dest_dir = Path(dest_dir)

    if not source.exists():
        raise FileNotFoundError(f"Source not found: {source}")

    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / source.name

    if dest.exists():
        if dest.is_dir():
            shutil.rmtree(dest)
        else:
            dest.unlink()

    shutil.move(str(source), str(dest))

    return dest


def s_common(parser):
    """subparser common arguments"""
    parser.add_argument("--outdir", help="Output diretory.", required=True)
    parser.add_argument("--sample", help="Sample name.", required=True)
    parser.add_argument("--thread", help=HELP_DICT["thread"], default=16)
    parser.add_argument("--debug", help=HELP_DICT["debug"], action="store_true")
    return parser


class Step:
    """
    Step class
    """

    def __init__(self, args, display_title=None):
        """
        display_title controls the section title in HTML report
        force thread <=20
        """
        sys.stderr.write(f"CeleScope version: {__version__} ")
        sys.stderr.write(f"Args: {args}\n")
        self.args = args
        self.outdir = args.outdir
        self.outs_dir = f"{args.outdir}/../outs"
        self.sample = args.sample
        self.assay = args.subparser_assay
        self.thread = int(args.thread)
        self.thread = min(self.thread, 20)
        self.debug = args.debug
        self.out_prefix = f"{self.outdir}/{self.sample}"
        self.display_title = display_title

        # metrics index
        self.metric_index = 0

        utils.check_mkdir(self.outdir)

        # set
        class_name = self.__class__.__name__
        if not display_title:
            self._display_title = class_name
        else:
            self._display_title = display_title
        self._step_name = class_name[0].lower() + class_name[1:]
        self.__slots = ["data", "metrics"]
        self._step_summary_name = f"{self._step_name}_summary"

        self.__metric_list = []
        self.__help_content = []
        self.__comments = []
        self._path_dict = {}
        for slot in self.__slots:
            self._path_dict[slot] = f"{self.outdir}/../.{slot}.json"

        self.__content_dict = {}
        self.old_step_dict = {}
        for slot, path in self._path_dict.items():
            if not os.path.exists(path):
                self.__content_dict[slot] = {}
                if slot == "data":
                    self.__content_dict[slot]["parameters"] = {}
            else:
                with open(path) as f:
                    try:
                        self.__content_dict[slot] = json.load(f)
                    except ValueError:
                        print(
                            f'WARNING: Decoding "{path}" as json has failed. Will create empty json file.'
                        )
                        self.__content_dict[slot] = {}
            # clear step_summary
            if self._step_summary_name in self.__content_dict[slot]:
                self.old_step_dict[slot] = self.__content_dict[slot][
                    self._step_summary_name
                ]
            self.__content_dict[slot][self._step_summary_name] = {}

        # jinja env
        self.env = Environment(
            loader=FileSystemLoader(os.path.dirname(__file__) + "/../templates/"),
            autoescape=select_autoescape(["html", "xml"]),
        )

        # out file
        self.__stat_file = f"{self.outdir}/stat.txt"
        self.report_html = f"{self.outdir}/../{self.sample}_report.html"

        # move file to outs
        self.outs = []

    def add_metric(
        self,
        name,
        value,
        total=None,
        help_info=None,
        display=None,
        show=True,
        print_log=True,
        value_type=None,
    ):
        """
        add metric to metric_list

        Args
            total: int or float, used to calculate fraction
            help_info: str, help info for metric in html report
            display: str, controls how to display the metric in HTML report.
            show: bool, whether to add to `.data.json` and `stat.txt`. `.data.json` is used for HTML report. `stat.txt` is used in house.
            print_log: bool, whether to print metric to stdout
            value_type: if value_type is "fraction", value = round(value * 100, 2), display = f'{value}%'
        """

        name = cap_str_except_preposition(name)
        if help_info:
            help_info = help_info[0].upper() + help_info[1:]
            if help_info[-1] != ".":
                help_info += "."
        if not display:
            if isinstance(value, numbers.Number):
                display = str(format(value, ","))
            else:
                display = value
        fraction = None
        if total:
            fraction = round(value / total * 100, 2)
            display += f"({fraction}%)"
        if value_type == "fraction":
            value = round(value * 100, 2)
            display = f"{value}%"
        self.__metric_list.append(
            {
                "name": name,
                "value": value,
                "total": total,
                "fraction": fraction,
                "display": display,
                "help_info": help_info,
                "show": show,
            }
        )

        if print_log:
            sys.stderr.write(f"{name}: {display}\n")

    def _write_stat(self):
        with open(self.__stat_file, "w") as writer:
            for metric in self.__metric_list:
                if metric["show"]:
                    name = metric["name"]
                    display = metric["display"]

                    line = f"{name}: {display}"
                    writer.write(line + "\n")

    def _dump_content(self):
        """dump content to json file"""
        for slot, path in self._path_dict.items():
            if self.__content_dict[slot]:
                with open(path, "w") as f:
                    json.dump(self.__content_dict[slot], f, indent=4)

    @utils.add_log
    def _render_html(self):
        template = self.env.get_template(f"html/{self.assay}/base.html")
        with io.open(self.report_html, "w", encoding="utf8") as f:
            html = template.render(self.__content_dict["data"])
            f.write(html)

    def _add_content_data(self):
        step_summary = {}
        step_summary["display_title"] = self._display_title
        metric_list = []
        comment_metric_list = []
        for metric in self.__metric_list:
            if metric["show"]:
                metric_list.append(metric)
            else:
                comment_metric_list.append(metric)
        step_summary["metric_list"] = metric_list
        step_summary["comment_metric_list"] = comment_metric_list
        step_summary["help_content"] = self.__help_content
        step_summary["comments"] = self.__comments
        self.__content_dict["data"][self._step_summary_name].update(step_summary)

    def _add_content_metric(self):
        metric_dict = dict()
        for metric in self.__metric_list:
            name = metric["name"]
            value = metric["value"]
            fraction = metric["fraction"]
            metric_dict[name] = value
            if fraction:
                metric_dict[f"{name} Fraction"] = fraction

        self.__content_dict["metrics"][self._step_summary_name].update(metric_dict)

    def add_data(self, **kwargs):
        """
        add data(other than metrics) to self.content_dict['data']
        for example: add plots and tables
        """
        for key, value in kwargs.items():
            self.__content_dict["data"][self._step_summary_name][key] = value

    def add_help_content(self, name, content):
        """
        add help info before metrics' help_info
        """
        self.__help_content.append({"name": name, "content": content})

    def add_comments(self, content):
        """
        add comment to self.__comments
        """
        self.__comments.append(content)

    @utils.add_log
    def get_slot_key(self, slot, step_name, key):
        """read slot from json file"""
        try:
            return self.__content_dict[slot][step_name + "_summary"][key]
        except KeyError:
            self.get_slot_key.logger.warning(
                f"{key} not found in {step_name}_summary.{slot}"
            )
            raise

    @utils.add_log
    def add_slot_step(self, slot, step_name, val):
        """add slot to json"""
        self.__content_dict[slot][step_name + "_summary"] = val

    @utils.add_log
    def get_slot_step(self, slot, step_name):
        """add slot to json"""
        return self.__content_dict[slot][step_name + "_summary"]

    @utils.add_log
    def add_table(self, title, table_id, df, help="", script=""):
        if not script:
            script = """
            <script>
            $(document).ready(function () {{
                var table = $('#{table_id}').DataTable({{
                    "order": [],
                    "scrollX": true,
                    "autoWidth": false,
                    dom: 'Bfrtip',
                    buttons: ['excel']
                }});
                table.draw();
            }});
            </script>
            """.format(table_id=table_id)

        self.__content_dict["data"][table_id] = {
            "title": title,
            "df": df.to_html(
                escape=False, index=False, table_id=table_id, justify="center"
            ),
            "help": help,
            "script": script,
        }

    def get_table_dict(self, title, table_id, df_table):
        """
        table_dict {title: '', table_id: '', df_table: pd.DataFrame}
        """
        table_dict = {}
        table_dict["title"] = title
        table_dict["table"] = df_table.to_html(
            escape=False, index=False, table_id=table_id, justify="center"
        )
        table_dict["id"] = table_id
        return table_dict

    def _move_files(self):
        for f in self.outs:
            try:
                safe_move(f, self.outs_dir)
            except FileNotFoundError as e:
                sys.stderr.write(f"WARNING: {e}\n")

    @utils.add_log
    def _add_parameters(self):
        if "parameters" not in self.__content_dict["data"]:
            sys.stderr.write("parameters not in {sample}/.data.json\n")
            self.__content_dict["data"]["parameters"] = {}
        for key, value in vars(self.args).items():
            if key not in {
                "func",
                "thread",
                "outdir",
                "sample",
                "subparser_assay",
                "debug",
            }:
                self.__content_dict["data"]["parameters"][key] = value

    @utils.add_log
    def _create_outs_link(self):
        if self.outs:
            os.makedirs(self.outs_dir, exist_ok=True)
        for f in self.outs:
            src = Path(f)
            dest = Path(self.outs_dir) / src.name
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            relative_src = os.path.relpath(src, start=dest.parent)
            os.symlink(relative_src, dest)

    @utils.add_log
    def _clean_up(self):
        self._add_content_data()
        self._add_content_metric()
        self._add_parameters()
        self._write_stat()
        self._dump_content()
        self._render_html()
        self._move_files()

    @utils.add_log
    def debug_subprocess_call(self, cmd):
        """
        debug subprocess call
        """
        self.debug_subprocess_call.logger.info(cmd)
        if cmd.find("2>&1") == -1:
            cmd += " 2>&1 "
        subprocess.check_call(cmd, shell=True)

    def get_metric_list(self):
        return self.__metric_list

    def set_metric_list(self, metric_list):
        self.__metric_list = metric_list

    @abc.abstractmethod
    def run(self):
        sys.exit("Please implement run() method.")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self._clean_up()
