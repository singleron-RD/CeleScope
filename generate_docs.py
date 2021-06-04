import argparse
import inspect
import os

import celescope.tools.utils as utils
from celescope.celescope import ArgFormatter
from celescope.__init__ import ASSAY_DICT

PRE_PROCESSING_STEPS = ('barcode', 'cutadapt')
DOCS_ROOT = 'docs'
PRE_PROCESSING_ROOT = 'pre_processing'


def generate_single_step_doc(assay, step):
    step_module = utils.find_step_module(assay, step)
    func_opts = getattr(step_module, f"get_opts_{step}")
    
    class_docs = get_class_docs(step_module)
    argument_docs = get_argument_docs(func_opts)

    if step in PRE_PROCESSING_STEPS:
        out_md = f'{DOCS_ROOT}/{PRE_PROCESSING_ROOT}/{step}.md'
    else:
        out_md = f'{DOCS_ROOT}/{assay}/{step}.md'

    with open(out_md, 'w') as out_file:
        out_file.write(class_docs)
        out_file.write(argument_docs)


def get_argument_docs(func_opts):
    argument_docs = ""
    parser = argparse.ArgumentParser(description='CeleScope',formatter_class=ArgFormatter)
    func_opts(parser, sub_program=True)
    for argument in parser._option_string_actions:
        if not argument in ['-h', '--help']:
            argument_docs += (f'`{argument}` {parser._option_string_actions[argument].help}\n\n')
    argument_docs = "## Arguments\n" + argument_docs
    return argument_docs


def get_class_docs(step_module):
    titles = ("Features", "Output")
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


@utils.add_log
def generate_all_docs():
    os.system(f'mkdir -p {DOCS_ROOT}/{PRE_PROCESSING_ROOT}')
    for assay in ASSAY_DICT:
        init_module = utils.find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__

        for step in __STEPS__:
            if not os.path.exists(assay):
                os.system(f'mkdir -p {DOCS_ROOT}/{assay}')
            generate_single_step_doc(assay, step)

if __name__ == "__main__":
    generate_all_docs()