from celescope.tools.make_ref import MakeRef_STAR


def mkref(args):
    genome_type = "virus"
    with MakeRef_STAR(genome_type, args) as runner:
        runner.run()


def get_opts_mkref(parser, sub_program):
    MakeRef_STAR.opts(parser, sub_program)
