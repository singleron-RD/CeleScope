import argparse
import importlib
from celescope.__init__ import __VERSION__, ASSAY_DICT


def main():
    """celescope cli
    """
    parser = argparse.ArgumentParser(description='CeleScope')
    parser.add_argument('-v', '--version', action='version', version=__VERSION__)
    subparsers = parser.add_subparsers()

    for assay in ASSAY_DICT:
        text = ASSAY_DICT[assay]
        subparser_1st = subparsers.add_parser(assay, description=text)
        # add 2ed subparser
        subparser_2ed = subparser_1st.add_subparsers()

        # import __STEPS__
        init_module = importlib.import_module(f"celescope.{assay}.__init__")
        __STEPS__ = init_module.__STEPS__

        for step in __STEPS__:
            # import function and opts
            try:
                step_module = importlib.import_module(f"celescope.{assay}.{step}")
            except ModuleNotFoundError as error:
                try:
                    step_module = importlib.import_module(f"celescope.tools.{step}")
                except ModuleNotFoundError as error:
                    module_path = init_module.IMPORT_DICT[step]
                    step_module = importlib.import_module(f"{module_path}.{step}")

            func = getattr(step_module, step)
            func_opts = getattr(step_module, f"get_opts_{step}")
            parser_step = subparser_2ed.add_parser(step)
            func_opts(parser_step, True)
            parser_step.set_defaults(func=func)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
