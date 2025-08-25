from celescope.ffpe.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_ffpe(Multi):
    pass


def main():
    multi = Multi_ffpe(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
