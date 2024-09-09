import argparse

# On some system, import scanpy will report "ImportError: dlopen: cannot load any more object with static TLS" when scanpy import sklearn
# as described in
# https://github.com/singleron-RD/CeleScope/issues/127
# https://github.com/scikit-learn/scikit-learn/issues/14485
# https://github.com/scverse/scanpy/issues/1121
# https://github.com/pytorch/pytorch/issues/2575
# After testing, the version of glibc and sklearn is not the reason
# Import scanpy before scipy.stats(used in tools/count.py) will do the trick.
# However, this approach is only tested on one system.
import scanpy  # noqa # pylint: disable=unused-import

from celescope.tools import utils
from celescope.__init__ import __VERSION__, ASSAY_LIST


class ArgFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


def main():
    """celescope cli"""
    parser = argparse.ArgumentParser(
        description="CeleScope", formatter_class=ArgFormatter
    )
    parser.add_argument("-v", "--version", action="version", version=__VERSION__)
    subparsers = parser.add_subparsers(dest="subparser_assay")

    for assay in ASSAY_LIST:
        subparser_1st = subparsers.add_parser(assay)

        # add 2ed subparser
        subparser_2nd = subparser_1st.add_subparsers()

        # import STEPS
        init_module = utils.find_assay_init(assay)
        STEPS = init_module.STEPS

        for step in STEPS:
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


if __name__ == "__main__":
    main()
