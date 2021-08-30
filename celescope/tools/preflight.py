import os
import subprocess
import scripts.release_local as local
import glob
import celescope.tools.utils as utils

'''
Functions for preflight checks 
'''

@utils.add_log
def check_file(mapfile):
    print("Checking mapfile format and raw fastq file path...")
    with open(mapfile) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line_split = line.split()
            library_id, library_path, sample_name = line_split[:3]
            if not library_path.startswith('/'):
                raise Exception("Specified %s file must be an absolute path: %s" % (library_id, library_path))
            if not os.path.exists(library_path):
                raise Exception("Specified %s file does not exist: %s" % (library_id, library_path))
            if not os.path.isdir(library_path):
                raise Exception("Specified %s file is a folder: %s" % (library_id, library_path))
            if not os.access(library_path, os.R_OK):
                raise Exception("Specified %s file is not readable: %s" % (library_id, library_path))

@utils.add_log
def check_soft(): 

    print("Checking required tools and packages for conda environment (%s) ..." % (os.path.basename(os.environ['CONDA_DEFAULT_ENV'])))
    # check Bin
    try:
        subprocess.check_call(["which", "STAR"])
    except subprocess.CalledProcessError:
        msg = "STAR not found on PATH. (Required for celescope)" 
        return msg
    try:
        subprocess.check_call(["which", "cutadapt"])
    except subprocess.CalledProcessError:
        msg = "cutadapt not found on PATH. (Required for celescope)" 
        return msg
    try:
        subprocess.check_call(["which", "featureCounts"])
    except subprocess.CalledProcessError:
        msg = "featureCounts not found on PATH. (Required for celescope)" 
        return msg
    try:
        subprocess.check_call(["which","picard"])
    except subprocess.CalledProcessError:
        msg = "picard not found on PATH. (Required for celescope)"
        return msg
    try:
        subprocess.check_call(["which","samtools"])
    except subprocess.CalledProcessError:
        msg = "samtools not found on PATH. (Required for celescope)"
        return msg
    try:
        subprocess.check_call(["which","mixcr"])
    except subprocess.CalledProcessError:
        msg = "mixcr not found on PATH. (Required for celescope)"
        return msg
    try:
        subprocess.check_call(["which","subread-align"])
    except subprocess.CalledProcessError:
        msg = "subread-align not found on PATH. (Required for celescope)"
        return msg
    try:
        subprocess.check_call(["which","subread-buildindex"])
    except subprocess.CalledProcessError:
        msg = "subread-buildindex not found on PATH. (Required for celescope)"
        return msg
    try:
        subprocess.check_call(["which","subread-fullscan"])
    except subprocess.CalledProcessError:
        msg = "subread-fullscan not found on PATH. (Required for celescope)"
        return msg
    try:
        subprocess.check_call(["which","gatk"])
    except subprocess.CalledProcessError:
        msg = "gatk not found on PATH. (Required for celescope)"
        return msg
    try:
        subprocess.check_call(["which","bcftools"])
    except subprocess.CalledProcessError:
        msg = "bcftools not found on PATH. (Required for celescope)"
        return msg

    # check R-package
    if not os.path.exists(os.environ['CONDA_DEFAULT_ENV']+'/lib/R/library/Seurat'):
        msg = "R package Seurat not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(os.environ['CONDA_DEFAULT_ENV']+'/lib/R/library/argparser'):
        msg = "R package argparser not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(os.environ['CONDA_DEFAULT_ENV']+'/lib/R/library/tidyverse'):
        msg = "R package tidyverse not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(os.environ['CONDA_DEFAULT_ENV']+'/lib/R/library/DropletUtils'):
        msg = "R package DropletUtils not found on PATH (Required for celescope)"
        raise Exception(msg)
    '''
    use for test error
    if not os.path.exists(os.environ['CONDA_DEFAULT_ENV']+'/lib/R/library/inferCNV'):
        msg = "R package inferCNV not found on PATH (Required for celescope)"
        raise Exception(msg)
    '''

    # check python package
    pypath = glob.glob(os.environ['CONDA_DEFAULT_ENV']+'/lib/*/site-packages/')[0]
    
    if not os.path.exists(pypath + 'pysam'):
        msg = "python package pysam not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'cutadapt'):
        msg = "python package cutadapt not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'scipy'):
        msg = "python package scipy not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'numpy'):
        msg = "python package numpy not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'xopen'):
        msg = "python package xopen not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'editdistance'):
        msg = "python package editdistance not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'sklearn'):
        msg = "python package sklearn not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'plotly'):
        msg = "python package plotly not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'plotnine'):
        msg = "python package plotnine not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'matplotlib'):
        msg = "python package matplotlib not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'Cython'):
        msg = "python package Cython not found on PATH (Required for celescope)"
        raise Exception(msg)
    if not os.path.exists(pypath + 'pytest'):
        msg = "python package pytest not found on PATH (Required for celescope)"
        raise Exception(msg)

@utils.add_log
def check_env():
    print("Checking conda environment (%s) ..." % (os.path.basename(os.environ['CONDA_DEFAULT_ENV'])))
    if local.ENV_NAME != os.path.basename(os.environ['CONDA_DEFAULT_ENV']):
        msg = "Wanring: Current conda environment's name doesn't match celescope version's name."
    else:
        msg = "You're in the right conda environment."
    msg += "\n" + "Make sure current environment has installed celescope and required packages successfully."
    print(msg)
