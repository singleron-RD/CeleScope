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

> [!NOTE] 
> Only **CeleScope v2.3.0 or higher** supports analysis of **GEXSCOPE-V3** data. For details, see [here](./chemistry.md).

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
|[bulk_vdj](assay/multi_bulk_vdj.md)|bulk VDJ|AccuraCode|
|[bulk_rna](assay/multi_bulk_rna.md)|bulk RNA|AccuraCode|
|[ffpe](assay/multi_ffpe.md)|single-cell FFPE||
|[space](assay/multi_space.md)|single-cell Spatial||


Click on the hyperlinks above to view the uasge for each assay. Run `multi_{assay} -h` in the command line to display available arguments. For example：
```
$ multi_rna -h

usage: rna multi-samples [-h] --mapfile MAPFILE [--mod {sjm,shell}] [--queue QUEUE] [--rm_files] [--steps_run STEPS_RUN]
...
```

To avoid dependency conflicts, certain workflows have been kept out of CeleScope and hosted in their own standalone repositories:

|data|kit|repo|
|---|------|--------------|
|single cell 3′ and 5′ RNA-seq|MobiuSCOPE|[celescope-mobiu](https://github.com/singleron-RD/celescope-mobiu)|
|single cell DNA + RNA|AccuraSCOPE|[AccuraSCOPE](https://github.com/singleron-RD/AccuraSCOPE)|
|single cell ATAC| |[CeleScope_ATAC](https://github.com/singleron-RD/CeleScope_ATAC)|

## [FAQ](./FAQ.md)

## [Chemistry](./chemistry.md)

## [Change log](./CHANGELOG.md)



 
