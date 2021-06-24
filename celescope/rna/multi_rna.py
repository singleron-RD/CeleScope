from celescope.rna.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_rna(Multi):
    pass


def main():
    multi = Multi_rna(__ASSAY__)
    multi.run()


if __name__ == '__main__':
    main()
