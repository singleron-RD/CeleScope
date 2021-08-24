import argparse
import inspect
import os

import celescope.tools.utils as utils
from celescope.celescope import ArgFormatter
from celescope.__init__ import ASSAY_DICT, RELEASED_ASSAYS, ROOT_PATH

PRE_PROCESSING_STEPS = ('sample', 'barcode', 'cutadapt')
DOCS_DIR = f'{ROOT_PATH}/../docs/'
TEMPLATE_DIR = 'docs_template/'
MANUAL = f'{DOCS_DIR}/manual.md'


def get_argument_docs_from_parser(parser):
    argument_docs = ""
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
        
        write_bool = False
        if doc:
            for title in titles:
                if title in doc:
                    write_bool = True

        if write_bool:
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
        self.steps = init_module.__STEPS__.copy()
        self.steps.append(f'multi_{assay}')

        self.out_md_dict = {}
        self.relative_md_path = {}
        self.release_bool = self.assay in RELEASED_ASSAYS

        assay_dir = f'docs/{assay}'
        if not os.path.exists(assay_dir):
            os.system(f'mkdir -p {assay_dir}')    

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


    def write_step_doc(self, step):
        """
        folder: docs/folder/*.md
        """
        step_module = utils.find_step_module(self.assay, step)
        folder = step_module.__name__.split('.')[1]
        self.out_md_dict[step] = f'docs/{folder}/{step}.md'
        self.relative_md_path[step] = f'{folder}/{step}.md'

        class_docs = get_class_docs(step_module)
        argument_docs = self.get_argument_docs(step, step_module)

        with open(self.out_md_dict[step], 'w') as out_file:
            out_file.write(class_docs)
            out_file.write(argument_docs)

    def run(self):
        if self.release_bool:
            with open(MANUAL, 'a') as writer:
                writer.write(f'## {ASSAY_DICT[self.assay]}\n')            

        for step in self.steps:
            self.write_step_doc(step)
            if self.release_bool:
                self.write_step_in_manual(step)


    def write_step_in_manual(self, step):
        """
        - [mkref](rna/mkref.md)
        """
        if not step in PRE_PROCESSING_STEPS:
            with open(MANUAL, 'a') as writer:
                writer.write(f'- [{step}]({self.relative_md_path[step]})\n')


def main():
    cmd = (
        f"rm -r {DOCS_DIR};"
        f"cp -r {TEMPLATE_DIR} {DOCS_DIR}"
    )
    os.system(cmd)

    for assay in ASSAY_DICT:
        docs_obj = Docs(assay)
        docs_obj.run()


if __name__ == "__main__":
    main()
