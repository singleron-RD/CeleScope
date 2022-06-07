import argparse
import inspect
import os
from collections import OrderedDict

from celescope.tools import utils
from celescope.celescope import ArgFormatter
from celescope.__init__ import ASSAY_LIST, RELEASED_ASSAYS

PRE_PROCESSING_STEPS = ('sample', 'barcode', 'cutadapt')
DOCS_DIR = f'docs/'
MANUAL = f'{DOCS_DIR}/manual.md'

def parse_docstring_to_dict(docstring):
    """
    Parse docstring(markown) to dict
    Return:
        dict {title: content}
    """

    doc_dict = {}
    title = None
    for line in docstring.split('\n'):
        if line.find('##') != -1:
            title = line.replace('##', '').strip().lower()
            doc_dict[title] = ""
        else:
            if title:
                doc_dict[title] += line + '\n'

    return doc_dict



def get_argument_docs_from_parser(parser):
    argument_docs = ""
    for argument in parser._option_string_actions:
        if not argument in ['-h', '--help']:
            help_msg = parser._option_string_actions[argument].help
            if help_msg:
                help_msg = help_msg.strip()
                # if end with code segment, do not add dot
                if help_msg[-1] != '.' and help_msg[-3:] != "```":
                    help_msg += '.'
            argument_docs += (f'`{argument}` {help_msg}\n\n')
    argument_docs = "## Arguments\n" + argument_docs
    return argument_docs


def get_class_docs(step_module):
    """
    Return:
        class_docs: str
        doc_dict: dict
    """
    # case ignored
    titles = ("features", "output", "usage")

    for child in inspect.getmembers(step_module, inspect.isclass):
        """Filter out class not defined in step_module"""
        class_obj = child[1]
        if class_obj.__module__ != step_module.__name__:
            continue
        doc = inspect.getdoc(class_obj)
        
        # only write class docs when the doc contains at least one title
        write_bool = False
        if doc:
            doc_dict = parse_docstring_to_dict(doc)
            for title in titles:
                if title in doc_dict:
                    write_bool = True
                    break

        if write_bool:
            class_docs = doc + '\n'
            return class_docs, doc_dict

    return "", {}


class Docs():
    def __init__(self, assay):
        self.assay = assay

        init_module = utils.find_assay_init(assay)
        self.steps = init_module.__STEPS__.copy()
        self.multi_step = f'multi_{assay}'

        self.out_md_dict = {}
        self.relative_md_path = {}
        self.release_bool = self.assay in RELEASED_ASSAYS

        # order: step
        self.output_file_doc_dict = OrderedDict()
        self.feature_doc_dict = OrderedDict()

        assay_dir = f'docs/{assay}'
        if not os.path.exists(assay_dir):
            os.system(f'mkdir -p {assay_dir}') 
        tools_dir = 'docs/tools'
        if not os.path.exists(tools_dir):
            os.system(f'mkdir {tools_dir}')

        # manual
        self.manual_lines = []


    @utils.add_log
    def get_argument_docs(self, step, step_module):
        self.get_argument_docs.logger.info(step)
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

    def get_multi_feature_doc(self):
        multi_feature_doc = "## Features\n"
        for step in self.feature_doc_dict:
            multi_feature_doc += f"### {step}\n"
            multi_feature_doc += self.feature_doc_dict[step]
            multi_feature_doc += '\n'

        return multi_feature_doc

    def get_multi_output_doc(self):
        multi_output_doc = "## Output files\n"
        for step in self.output_file_doc_dict:
            multi_output_doc += f"### {step}\n"
            multi_output_doc += self.output_file_doc_dict[step]
            multi_output_doc += '\n'

        return multi_output_doc

    def write_step_doc(self, step):
        """
        folder: docs/folder/*.md
        """
        step_module = utils.find_step_module(self.assay, step)
        folder = step_module.__name__.split('.')[1]
        self.out_md_dict[step] = f'docs/{folder}/{step}.md'
        self.relative_md_path[step] = f'{folder}/{step}.md'

        class_docs, doc_dict = get_class_docs(step_module)
        if 'output' in doc_dict:
            self.output_file_doc_dict[step] = doc_dict['output']
        if 'features' in doc_dict:
            self.feature_doc_dict[step] = doc_dict['features']
        argument_docs = self.get_argument_docs(step, step_module)

        with open(self.out_md_dict[step], 'w') as out_file:
            out_file.write(class_docs)
            if step.startswith("multi"):
                out_file.write(self.get_multi_feature_doc())
                out_file.write(self.get_multi_output_doc())
            out_file.write(argument_docs)

    def run(self):
        if self.release_bool:
            assay_text = utils.get_assay_text(self.assay)
            self.manual_lines.append(f'## {assay_text}\n')       

        for step in self.steps:
            self.write_step_doc(step)
            if self.release_bool:
                self.add_step_in_manual(step)

        # write multi last because multi need to collect output files from each step
        self.write_step_doc(self.multi_step)
        if self.release_bool:
            self.add_step_in_manual(self.multi_step)

        self.write_manual()


    def add_step_in_manual(self, step):
        """
        - [mkref](rna/mkref.md)
        """
        if not step in PRE_PROCESSING_STEPS:
            self.manual_lines.append(f'- [{step}]({self.relative_md_path[step]})\n')

    def write_manual(self):
        with open(MANUAL, 'a') as writer:
            for line in self.manual_lines:
                writer.write(line)


def remove_old_assays():
    contents = []
    with open(MANUAL, 'r') as reader:
        for line in reader:
            contents.append(line)
            if line.find("tools/cutadapt.md") != -1:
                break
    
    with open(MANUAL, 'w') as writer:
        for line in contents:
            writer.write(line)


def main():

    remove_old_assays()
    
    for assay in ASSAY_LIST:
        docs_obj = Docs(assay)
        docs_obj.run()


if __name__ == "__main__":
    main()
