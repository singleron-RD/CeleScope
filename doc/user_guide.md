# User guide

## Installation
1. Clone repo
```
git clone https://github.com/singleron-RD/CeleScope.git
```

2. Create conda environment and install conda packages. 
It is recommended to use [mamba](https://github.com/mamba-org/mamba) (which is a faster replacement for Conda):
```
conda install mamba
cd CeleScope
mamba create -n celescope -y --file conda_pkgs.txt
```


3. Install celescope

Make sure you have activated the `celescope` conda environment before running `pip install celescope`. 
```
conda activate celescope
pip install celescope
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
|[sweetseq](assay/multi_citeseq.md)|single-cell glycosylation|ProMoSCOPE<sup>TM</sup>|
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


 
