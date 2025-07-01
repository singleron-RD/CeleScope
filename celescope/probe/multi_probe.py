from celescope.probe.__init__ import __ASSAY__
from celescope.tools.multi import Multi


class Multi_probe(Multi):
    """
    probe and primer detection
    """


def main():
    multi = Multi_probe(__ASSAY__)
    multi.run()


if __name__ == "__main__":
    main()
