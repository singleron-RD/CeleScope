
# CeleScope
CeleScope is a collection of bioinfomatics analysis pipelines developed at Singleron to process single cell sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis as well as a summary of QC criteria.

Detailed docs can be found in [manual](./docs/manual.md).

## Hardware/Software Requirements

- minimum 32GB RAM(to run STAR aligner)
- conda
- git

## Installation

1. Clone repo
```
git clone https://github.com/singleron-RD/CeleScope.git
```

2. Create conda environment and install conda packages
```
cd CeleScope
conda create -n celescope
conda activate celescope
conda install -y --file conda_pkgs.txt
```

3. Install celescope

Make sure you have activated the `celescope` conda environment before running `pip install celescope`. 
```
pip install celescope
```

## [Quick start](./docs/quick_start.md)


