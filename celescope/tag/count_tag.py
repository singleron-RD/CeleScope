import argparse
from .Count_tag import Count_tag
from celescope.tools.Step import Step, s_common


def get_opts_count_tag(parser, sub_program):
    parser.add_argument(
        "--UMI_min",
        help="cells have tag_UMI>=UMI_min are considered as valid cell",
        default="auto")
    parser.add_argument("--dim", help="tag dimension", default=1)
    parser.add_argument(
        "--SNR_min",
        help="minimum signal to noise ratio",
        default="auto")
    parser.add_argument(
        "--combine_cluster",
        help="conbine cluster tsv file",
        default=None)
    parser.add_argument("--coefficient", "-c", help="SNR coefficient", default=0.1)
    if sub_program:
        parser = s_common(parser)
        parser.add_argument("--read_count_file", help="tag read count file")
        parser.add_argument("--match_dir", help="matched scRNA-Seq CeleScope directory path", required=True)


def count_tag(args):

    step_name = "count_tag"
    step = Step(args, step_name)

    count_tag_object = Count_tag(
        args.sample,
        args.outdir,
        args.assay,
        args.read_count_file,
        args.match_dir,
        args.UMI_min,
        args.SNR_min,
        args.combine_cluster,
        args.dim,
        args.coefficient
    )
    count_tag_object.run()
    count_tag_object.stat()

    step.clean_up()