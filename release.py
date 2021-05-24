import subprocess

from celescope.__init__ import __version__

ENV_NAME = f'celescope{__version__}'

def create_conda():
    cmd = f"""
    set -e
    conda create -n {ENV_NAME}
    source activate {ENV_NAME}
    conda install --file conda_pkgs.txt --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing

    pip install -i https://pypi.tuna.tsinghua.edu.cn/simple celescope
    python setup.py install
    """
    subprocess.check_call(cmd, shell=True)


def lint_code():
    cmd = """
    set -e
    celescope -h
    # lint
        # W1618 (no-absolute-import)
        # E1101 (no-member)
        # W1633 (round-builtin)
        # W1619 (old-division)
        # W0105 (String statement has no effect)
        # W0511 TODO!
        # E1130 bad operand type for unary ~: _isnan (invalid-unary-operand-type)
    pylint --disable=all --enable=E,W --disable=W1618,E1101,W1633,W1619,W0105,W0511,E1130 celescope
    """
    subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    lint_code()
    create_conda()
