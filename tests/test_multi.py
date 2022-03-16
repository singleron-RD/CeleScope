"""
Integration tests
"""

import os
import subprocess
from concurrent import futures

from celescope.tools import utils

ASSAYS = [
    'fusion',
    'vdj',
    'tag',
    'capture_virus',
    'snp',
    'rna',
    'dynaseq',
]


def run_single(assay, test_dir):
    """
    Returns:
        string indicates complete status
    """
    os.chdir(os.path.join(test_dir, assay))
    print("*" * 20 + "running " + assay + "*" * 20)
    subprocess.check_call('sh run_shell.sh', shell=True)
    subprocess.check_call('sh sjm.sh', shell=True)
    try:
        subprocess.check_call('sh ./shell/test1.sh', shell=True)
    except subprocess.CalledProcessError:
        return f"{assay} failed"
    print("*" * 20 + "success " + assay + "*" * 20)
    return f"{assay} success."


@utils.add_log
def test_mutiple(assays, test_dir):
    """
    Run all
    >>> pytest -s celescope/tests/test_multi.py --test_dir {some_dir}
    Run some tests
    >>> pytest -s celescope/tests/test_multi.py --test_dir {some_dir} --assays tag,fusion
    """

    if not assays:
        assays = ASSAYS
    else:
        assays = assays.split(',')
    print("assays to run: ", assays)
    thread = len(assays)
    executor = futures.ProcessPoolExecutor(max_workers=thread)
    results = executor.map(run_single, assays, [test_dir] * len(assays))
    res_list = []
    for result in results:
        res_list.append(result)
    for result in res_list:
        print(result)
    assert not any((string.find("failed") != -1 for string in res_list))
