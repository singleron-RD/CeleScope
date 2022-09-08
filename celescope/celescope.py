import argparse
import os
import celescope
import subprocess

from celescope.tools import utils
from celescope.__init__ import __VERSION__, ASSAY_LIST

TOOLS_DIR = os.path.dirname(celescope.tools.__file__)
SINGLE_MODULE = ['mkgtf']
ASSAY_LIST += SINGLE_MODULE

class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main():
    """celescope cli
    """
    parser = argparse.ArgumentParser(description='CeleScope', formatter_class=ArgFormatter)
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers(dest='subparser_assay')

    for assay in ASSAY_LIST:
        text = utils.get_assay_text(assay)
        subparser_1st = subparsers.add_parser(assay, description=text)

        if assay == 'mkgtf':
            def func(args):
                subprocess.check_call(f'python {TOOLS_DIR}/mkgtf.py {args.gtf} {args.filtered_gtf}', shell=True)
            
            subparser_1st.add_argument('gtf', help='raw gtf file')
            subparser_1st.add_argument('filtered_gtf', help='filtered gtf file')
            subparser_1st.set_defaults(func=func)
            continue

        # add 2ed subparser
        subparser_2nd = subparser_1st.add_subparsers()

        # import __STEPS__
        init_module = utils.find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__

        for step in __STEPS__:
            # import function and opts
            step_module = utils.find_step_module(assay, step)
            func = getattr(step_module, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            parser_step = subparser_2nd.add_parser(step, formatter_class=ArgFormatter)
            func_opts(parser_step, sub_program=True)
            parser_step.set_defaults(func=func)

    args = parser.parse_args()
    if len(args.__dict__) <= 1:
        # No arguments or subcommands were given.
        parser.print_help()
        parser.exit()
    else:
        args.func(args)


if __name__ == '__main__':
    main()
