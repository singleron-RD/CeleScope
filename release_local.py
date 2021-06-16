import subprocess

from celescope.__init__ import __version__
from celescope.tools.utils import add_log

ENV_NAME = f'celescope{__version__}'
CONDA_ROOT = '/SGRNJ/Public/Software/conda_env/'

@add_log
def create_conda():
    cmd = f"""
    set -eo pipefail
    conda create -n {ENV_NAME}
    source activate {ENV_NAME}
    conda install --file conda_pkgs.txt --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing

    pip install -i https://pypi.mirrors.ustc.edu.cn/simple/ celescope
    python setup.py install
    ln -s /SGRNJ/Database/script/soft/gatk-4.1.8.1/gatk {CONDA_ROOT}/{ENV_NAME}/bin/gatk
    """
    print(cmd)
    subprocess.check_call(cmd, shell=True)


@add_log
def lint_code():
    cmd = """
    set -eo pipefail
    celescope -h
    pip install -i https://pypi.mirrors.ustc.edu.cn/simple/ pylint
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
    print(cmd)
    subprocess.check_call(cmd, shell=True)


@add_log
def zip_wdl():
    cmd = "cd wdl/ && zip -r wdl.zip ./*"
    print(cmd)
    subprocess.check_call(cmd, shell=True)

@add_log
def test_wdl():
    cmd = (
        "cd /SGRNJ03/randd/user/zhouyiqi/temp/wdl; "
        "sh /SGRNJ/Database/script/pipe/develop/dev_CeleScope/wdl/rna/local/run.sh "
    )
    print(cmd)
    subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    lint_code()
    zip_wdl()
    test_wdl()
    create_conda()
