import argparse
import inspect
import os
import importlib
from collections import defaultdict

import celescope.tools.utils as utils
from celescope.celescope import ArgFormatter
from celescope.__init__ import ASSAY_DICT, RELEASED_ASSAYS

PRE_PROCESSING_STEPS = ('sample', 'barcode', 'cutadapt')
DOCS_DIR = 'docs/'
TEMPLATE_DIR = 'docs_template/'
MANUAL_MD = f'{DOCS_DIR}/manual.md'
MANUAL_TEMPLATE = f'{DOCS_DIR}/manual_template.md'


def get_argument_docs_from_parser(parser):
    for argument in parser._option_string_actions:
        if not argument in ['-h', '--help']:
            help_msg = parser._option_string_actions[argument].help
            if help_msg:
                help_msg = help_msg.strip()
            argument_docs += (f'`{argument}` {help_msg}\n\n')
    argument_docs = "## Arguments\n" + argument_docs
    return argument_docs


def get_class_docs(step_module):
    titles = ("Features", "Output", "Usage")
    class_docs = ""
    for child in inspect.getmembers(step_module, inspect.isclass):
        """Filter out class not defined in step_module"""
        class_obj = child[1]
        if class_obj.__module__ != step_module.__name__:
            continue
        doc = inspect.getdoc(class_obj)
        if doc and "Features" in doc:
            for line in doc.split('\n'):
                for title in titles:
                    if line.find(title) != -1:
                        class_docs += "## " + line + '\n'
                        break
                else:
                    class_docs += line + '\n'
    class_docs += '\n\n'
    return class_docs


class Docs():
    def __init__(self, assay):
        self.assay = assay

        init_module = utils.find_assay_init(assay)
        self.steps = init_module.__STEPS__
        self.steps.append(f'multi_{assay}')
        folder = f'{DOCS_DIR}/{assay}/'

        self.out_md_dict = {}
        self.relative_md_path = {}
        for step in self.steps:
            self.out_md_dict[step] = f'{folder}/{step}.md'
            self.relative_md_path[step] = f'{assay}/{step}.md'

        if not os.path.exists(folder):
            os.system(f'mkdir -p {folder}')    

    def get_argument_docs(self, step, step_module):
        if step.startswith("multi"):
            multi_class = getattr(step_module, f'Multi_{self.assay}')
            multi_obj = multi_class(self.assay)
            argument_docs = get_argument_docs_from_parser(multi_obj.parser)
        else:
            parser = argparse.ArgumentParser(description='CeleScope', formatter_class=ArgFormatter)
            func_opts = getattr(step_module, f"get_opts_{step}")
            func_opts(parser, sub_program=True)
            argument_docs = get_argument_docs_from_parser(parser)
        return argument_docs   


    def write_step_doc(self, step):
        step_module = utils.find_step_module(self.assay, step)
        class_docs = get_class_docs(step_module)
        argument_docs = self.get_argument_docs(step, step_module)

        with open(self.out_md_dict[step], 'w') as out_file:
            out_file.write(class_docs)
            out_file.write(argument_docs)

def write_step_in_manual(md_path, step, manual_handle):
    """
    - [mkref](rna/mkref.md)
    """
    if not step in PRE_PROCESSING_STEPS:
        manual_handle.write(f'- [{step}]({md_path})\n')
    

"""
@utils.add_log
def generate_all_docs():
    md_path_dict = defaultdict(dict)

    for assay in ASSAY_DICT:
        init_module = utils.find_assay_init(assay)
        steps = init_module.__STEPS__
        generate_all_docs.logger.info(f"Writing docs {assay} ")

        steps.append(f'multi_{assay}')
        for step in steps:
            generate_all_docs.logger.info(f"Writing doc {assay}.{step}")
            md_path = generate_single_step_doc(assay, step)
            md_path_dict[assay][step] = md_path
    return md_path_dict
"""


@utils.add_log
def write_manual(md_path_dict):
    with open(MANUAL_MD, 'w') as manual_handle:
        with open(MANUAL_TEMPLATE, 'r') as manual_template:
            manual_handle.write(manual_template.read())
        for assay in RELEASED_ASSAYS:
            init_module = utils.find_assay_init(assay)
            steps = init_module.__STEPS__
            manual_handle.write(f'## {ASSAY_DICT[assay]}\n')
            for step in steps:
                write_manual.logger.info(f"Writing manual {assay}.{step}")
                md_path = md_path_dict[assay][step]
                write_step_in_manual(md_path, step, manual_handle)



if __name__ == "__main__":
    cmd = f"cp -r {TEMPLATE_DIR} {DOCS_DIR}"
    os.system(cmd)
    