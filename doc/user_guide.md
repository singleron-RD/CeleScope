# User guide

## Installation

### Create conda environment and install conda packages. 
First, you need to get the txt file from the github repository containing the name of the conda package you need to install. You can download it directly from the github repository interface, or use a download link.
The following command will download the conda_pkgs.txt required for the latest version.
```
wget https://raw.githubusercontent.com/singleron-RD/CeleScope/master/conda_pkgs.txt
```
Use the mirror if you have trouble connecting to github:
```
wget https://gitee.com/singleron-rd/celescope/raw/master/conda_pkgs.txt
```
If you need to download an earlier version of `conda_pkgs.txt`, just change `master` in the link to that version
```
wget https://raw.githubusercontent.com/singleron-RD/CeleScope/v1.15.0/conda_pkgs.txt
```

Then, start creating the conda environment and install the conda package.It is recommended to use [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (which is a faster replacement for Conda) to install conda packages.
The following command will create a conda environment named `celescope` and install the dependency packages contained in conda_pkgs.txt.
```
mamba create -n celescope -y --file conda_pkgs.txt
```

### Install celescope

Make sure you have activated the conda environment before running `pip install celescope`. The following command will install the latest version of celescope.
```
conda activate celescope
pip install celescope
```

If you want to install a previous version instead of the latest version, for example v1.15.0
```
pip install celescope==1.15.0
```

## Usage

CeleScope contains interfaces `multi_{assay}` to generate pipeline scripts for all assays. Assays can be one of:
|assay|data|kit|
|---|------|--------------|
|[rna](./assay/multi_rna.md)|single-cell RNA-Seq|GEXSCOPE<sup>TM</sup>|
|[dynaseq](./assay/multi_dynaseq.md)|single-cell dynamic RNA-Seq|DynaSCOPE<sup>TM</sup>|
|[tag](./assay/multi_tag.md)|single-cell sample multiplexing|CLindex<sup>TM</sup>|
|[vdj](./assay/multi_vdj.md)|single-cell VDJ|GEXSCOPE<sup>TM</sup> IR|
|[flv_trust4](./assay/multi_flv_trust4.md)|single-cell full length VDJ|sCircle<sup>TM</sup>|
|[capture_virus](./assay/multi_capture_virus.md)|single-cell virus|FocuSCOPE<sup>TM</sup> mRNA × EBV|
|[snp](./assay/multi_snp.md)|single-cell variant|FocuSCOPE<sup>TM</sup>|
|[fusion](./assay/multi_fusion.md)|single-cell fusion|FocuSCOPE<sup>TM</sup>|
|[sweetseq](assay/multi_sweetseq.md)|single-cell glycosylation|ProMoSCOPE<sup>TM</sup>|
|[citeseq](assay/multi_citeseq.md)|single-cell CITE-Seq|NA|
|[bulk_vdj](assay/multi_bulk_vdj.md)|bulk_vdj|NA|


Click on the hyperlinks above to view the uasge for each assay. Run `multi_{assay} -h` in the command line to display available arguments. For example：
```
$ multi_rna -h

usage: rna multi-samples [-h] --mapfile MAPFILE [--mod {sjm,shell}] [--queue QUEUE] [--rm_files] [--steps_run STEPS_RUN]
...
```

## [Test scripts and data](https://github.com/singleron-RD/celescope_test_script)

## [Change log](./CHANGELOG.md)


 
