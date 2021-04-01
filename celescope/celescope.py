import argparse
import importlib
from celescope.tools.utils import *
from celescope.__init__ import __VERSION__, ASSAY_DICT


def main():
    """celescope cli
    """
    parser = argparse.ArgumentParser(description='CeleScope', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers()

    for assay in ASSAY_DICT:
        text = ASSAY_DICT[assay]
        subparser_1st = subparsers.add_parser(assay, description=text)
        # add 2ed subparser
        subparser_2ed = subparser_1st.add_subparsers()

        # import __STEPS__
        init_module = find_assay_init(assay)
        __STEPS__ = init_module.__STEPS__

        for step in __STEPS__:
            # import function and opts
            step_module = find_step_module(assay, step)
            func = getattr(step_module, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            parser_step = subparser_2ed.add_parser(step, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            func_opts(parser_step, sub_program=True)
            parser_step.set_defaults(func=func)


    # vdj10X
    assay = 'vdj10X'
    text = ASSAY_DICT[assay]
    subparsers_assay = subparsers.add_parser(
        assay, help=text, description=text)
    subparsers_assay_sub = subparsers_assay.add_subparsers()

    parser_tmp = subparsers_assay_sub.add_parser('sample')
    get_opts_sample(parser_tmp, True)
    parser_tmp.set_defaults(func=sample_info)

    from celescope.vdj10X.convert import convert, get_opts_convert
    parser_tmp = subparsers_assay_sub.add_parser('convert')
    get_opts_convert(parser_tmp, True)
    parser_tmp.set_defaults(func=convert)

    from celescope.vdj10X.vdj_10X import vdj_10X, get_opts_vdj_10X
    parser_tmp = subparsers_assay_sub.add_parser('vdj_10X')
    get_opts_vdj_10X(parser_tmp, True)
    parser_tmp.set_defaults(func=vdj_10X)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
