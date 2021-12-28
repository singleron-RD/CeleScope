import argparse

import celescope.tools.utils as utils
from celescope.__init__ import __VERSION__, ASSAY_LIST


class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def main():
    """celescope cli
    """
    parser = argparse.ArgumentParser(description='CeleScope', formatter_class=ArgFormatter)
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers()

    for assay in ASSAY_LIST:
        text = utils.get_assay_text(assay)
        subparser_1st = subparsers.add_parser(assay, description=text)
        # add 2ed subparser
        subparser_2ed = subparser_1st.add_subparsers()

        # import __STEPS__
        init_module = utils.find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__

        for step in __STEPS__:
            # import function and opts
            step_module = utils.find_step_module(assay, step)
            func = getattr(step_module, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            parser_step = subparser_2ed.add_parser(step, formatter_class=ArgFormatter)
            func_opts(parser_step, sub_program=True)
            parser_step.set_defaults(func=func)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
