import scripts.release_local as local
import celescope.tools.utils as utils

import glob
import os
import subprocess
import re
import pandas as pd


'''
Functions for preflight checks 

celescope required package vision
'''
ref_vers ={'python': '3.6.7', 
            'star': '2.6.1b', 
            'cutadapt': '1.17',
            'samtools': '1.9', 
            'picard': '2.18.17', 
            'mixcr': '3.0.3', 
            'subread': '2.0.1', 
            'pysam': '0.16.0.1', 
            'scipy': '1.5.0', 
            'numpy': '1.19.5', 
            'xopen': '0.5.0', # >=
            'sklearn': '0.0', 
            'plotly': '4.14.3', 
            'matplotlib': '3.3.0', 
            'plotnine': '0.8.0', 
            'editdistance': '0.5.3', #>=
            'Seurat': '4.0.1',
            'jinja2': '2.10'} #>=


@utils.add_log
def check_soft(): 
    print("Checking required tools and packages for environment (%s) ..." % (os.path.basename(os.environ['CONDA_DEFAULT_ENV'])))
    
    # check Bins
    Bins = ['STAR','cutadapt','featureCounts','picard','samtools','mixcr',
           'subread-align','subread-buildindex','subread-fullscan','gatk','bcftools']
    for tool in Bins:
        try:
            subprocess.check_output(['which',tool])
        except subprocess.CalledProcessError:
            msg = tool + " not found on PATH. (Required for celescope)"
            return msg

    # check R-packages
    R_pack = ['Seurat','argparser','tidyverse','DropletUtils']
    for pack in R_pack:
        path = os.environ['CONDA_DEFAULT_ENV'] + '/lib/R/library/' + pack
        if not os.path.exists(path):
            msg = "R package " + pack + " not found on PATH (Required for celescope)"
            raise Exception(msg)

    # check python packages
    pypack = ['pysam','cutadapt','scipy','numpy','xopen',
              'editdistance','sklearn','plotly','plotnine','matplotlib','Cython','pytest']
    
    pypath = glob.glob(os.environ['CONDA_DEFAULT_ENV']+'/lib/*/site-packages/')[0]
    
    for pack in pypack:
        path = pypath + pack
        if not os.path.exists(path):
            msg = "python package " + pack + " not found on PATH (Required for celescope)"
            raise Exception(msg)


@utils.add_log
def check_env():
    print("Checking Environment (%s) ..." % (os.path.basename(os.environ['CONDA_DEFAULT_ENV'])))
    
    if local.ENV_NAME != os.path.basename(os.environ['CONDA_DEFAULT_ENV']):
        msg = "Warning: Current conda environment's name doesn't match celescope version's name."
    else:
        msg = "You're in the right conda environment."
    msg += "\n" + "         Make sure celescope and required packages have been installed in current environment successfully."
    print(msg)


def get_version():
    '''
    get version cmd
    '''
    PACKAGE_VERSION_CMDS = []
    
    # Bin
    cmd = 'conda list'
    res = subprocess.check_output(cmd,shell=True)
    res = res.decode('utf8').strip().split()
    Bins = ['python','star','cutadapt','samtools','picard','mixcr','subread','bcftools','featureCounts']
    for Bin in Bins:
        if Bin in res:
            PACKAGE_VERSION_CMDS += [{'name':Bin, 'version':res[res.index(Bin)+1]}]
        elif Bin == 'featureCounts':
            cmd = 'featureCounts 2>&1| grep "^ *Version"' 
            version = subprocess.check_output(cmd, shell=True)
            version = version.decode('utf8').strip().split()[1]
            PACKAGE_VERSION_CMDS += [{'name':Bin, 'version':version}]

    # py-package
    pypacks = ['pysam','cutadapt','scipy','numpy','xopen','jinja2','sklearn','plotly','matplotlib','Cython','pytest','plotnine','editdistance']
    for pack in pypacks:        
        cmd = 'python -c "import %s ;print(%s.__version__)" ' % (pack,pack)
        PACKAGE_VERSION_CMDS += [{'name':pack,'cmd':cmd}]
    
    # R package 
    R_packs = ['Seurat','argparser','tidyverse','DropletUtils']
    for pack in R_packs:
        version = os.environ['CONDA_DEFAULT_ENV'] + '/lib/R/library/' + pack + '/DESCRIPTION'
        cmd = 'cat ' + version + '|grep "^ *Version"'
        PACKAGE_VERSION_CMDS += [{'name':pack,'cmd':cmd}]
    
    return PACKAGE_VERSION_CMDS


@utils.add_log
def record_package_versions(outdir,PACKAGE_VERSION_CMDS):
    print("Checking required tools and packages VERSION for environment (%s) ..." % (os.path.basename(os.environ['CONDA_DEFAULT_ENV'])))
    head = []
    data = []
    for package in PACKAGE_VERSION_CMDS:
        name = package['name']
        if 'cmd' in package.keys():
            cmd = package['cmd']
            version = "not found"
            
            try:
                version = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
                version = version.decode('utf8').strip()
                if version.startswith('Version'):
                    version = version.split()[1]

            except subprocess.CalledProcessError:
                #  has no attribute '__version__'
                print('package %s has no attribute __version__, finding %s in packages directory...'%(name,name))
                pypath = glob.glob(os.environ['CONDA_DEFAULT_ENV']+'/lib/*/site-packages/')[0]
                files = os.listdir(pypath)
                for file in files:
                    if name in file:
                        version = re.findall(r'\d+',file)
                        version = '.'.join(version)
                    else:
                        continue
            head.append(name)
            data.append(version)
        else:
            head.append(name)
            data.append(package['version'])

    ver_dict = dict(zip(head,data))

    # Write to Version.txt file
    p_list = []
    for i in range(len(head)):
        p_list.append({"name":head[i], "version":data[i]})
    p_df = pd.DataFrame(p_list)
    p_df.to_csv(f'{outdir}/Version.txt', sep='\t', index=False, header=None)
    return ver_dict


def compare_version(ver_dict,ref_dict):
    '''
    Compare package version to celescope required package version
    '''
    msg = "package %s version is %s, celescope required %s version is >= %s" 
    for key,val in ver_dict.items():
        if key == 'xopen':
            if val != ref_dict[key]:
                print(msg%(key,ver_dict[key],key,ref_dict[key]))
        elif key == 'jinja2':
            if val != ref_dict[key]:
                print(msg%(key,ver_dict[key],key,ref_dict[key]))
        elif key == 'editdistance':
            if val != ref_dict[key]:
                print(msg%(key,ver_dict[key],key,ref_dict[key]))
        elif key in ref_dict:
            if ver_dict[key] != ref_dict[key]:
                print("Warning: package %s version is %s, which doesn't match celescope required version: %s" 
                      %(key,ver_dict[key],ref_dict[key]))


def run_preflight(outdir):
    check_env()
    check_soft()
    ver_dict = record_package_versions(outdir,PACKAGE_VERSION_CMDS=get_version())
    compare_version(ver_dict,ref_dict=ref_vers)
