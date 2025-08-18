# User guide

## Requirements

- 64 bit Linux
- Minimum 32GB RAM to run [STAR](https://github.com/alexdobin/STAR) with human/mouse genome


## Installation

### Docker image

https://quay.io/repository/singleron-rd/celescope?tab=tags

If you are unable to use Docker, you can install it using mamba/conda and pip as follows.

### Create conda environment and install conda packages. 
First, you need to get the txt file from the github repository containing the name of the conda package you need to install. You can download it directly from the github repository interface, or use a download link.
The following command will download the conda_pkgs.txt required for the latest version.
```
wget https://raw.githubusercontent.com/singleron-RD/CeleScope/master/conda_pkgs.txt
```

Note: Accessing raw files from GitHub can be unstable within mainland China. A Gitee mirror is provided for more reliable downloading:
```
wget https://gitee.com/singleron-rd/celescope/raw/master/conda_pkgs.txt
```


Then, start creating the conda environment and install the conda package.It is recommended to use [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (which is a faster replacement for Conda) to install conda packages.
The following command will create a conda environment named `celescope` and install the dependency packages contained in conda_pkgs.txt.
```
mamba create -n celescope -y --file conda_pkgs.txt
```

### Install celescope

Make sure you have activated the conda environment before running `pip install celescope`. The following command will install the latest version of celescope.
```
mamba activate celescope
pip install celescope
```

## Usage

CeleScope contains interfaces `multi_{assay}` to generate pipeline scripts for all assays. Assays can be one of:
|assay|data|kit|
|---|------|--------------|
|[rna](./assay/multi_rna.md)|single-cell RNA-Seq|GEXSCOPE|
|[dynaseq](./assay/multi_dynaseq.md)|single-cell dynamic RNA-Seq|DynaSCOPE|
|[tag](./assay/multi_tag.md)|single-cell sample multiplexing|CLindex|
|[vdj](./assay/multi_vdj.md)|single-cell VDJ|GEXSCOPE IR|
|[flv_trust4](./assay/multi_flv_trust4.md)|single-cell full length VDJ|sCircle|
|[capture_virus](./assay/multi_capture_virus.md)|single-cell virus|FocuSCOPE mRNA × EBV|
|[snp](./assay/multi_snp.md)|single-cell variant|FocuSCOPE|
|[fusion](./assay/multi_fusion.md)|single-cell fusion|FocuSCOPE|
|[sweetseq](assay/multi_sweetseq.md)|single-cell glycosylation|ProMoSCOPE|
|[citeseq](assay/multi_citeseq.md)|single-cell CITE-Seq|NA|
|[bulk_vdj](assay/multi_bulk_vdj.md)|bulk_vdj|AccuraCode|
|[bulk_rna](assay/multi_bulk_rna.md)|bulk_rna|AccuraCode|


Click on the hyperlinks above to view the uasge for each assay. Run `multi_{assay} -h` in the command line to display available arguments. For example：
```
$ multi_rna -h

usage: rna multi-samples [-h] --mapfile MAPFILE [--mod {sjm,shell}] [--queue QUEUE] [--rm_files] [--steps_run STEPS_RUN]
...
```


## [Test scripts and data](https://github.com/singleron-RD/celescope_test_script)

## [Change log](./CHANGELOG.md)

## Frequently asked questions

### What if CeleScope reports an error?

When CeleScope reports an error, using GitHub Issues can be an effective way to solve or report the error.

Before creating a new issue, search the existing issues with keywords(e.g., error message) to see if someone else has already reported the same problem. 

### How to install a previous version instead of the latest version?

1. Download an earlier version of `conda_pkgs.txt`: just change `master` to the version number(e.g. v.15.0) in the link to that version 
```
wget https://raw.githubusercontent.com/singleron-RD/CeleScope/v1.15.0/conda_pkgs.txt
# create the env
mamba create -n celescope1.15.0 -y --file conda_pkgs.txt
```

2. Install celescope
```
mamba activate celescope1.15.0
pip install celescope==1.15.0
```



 
