# Hardware/Software Requirements

- minimum 32GB RAM(to run STAR aligner)
- conda
- git

# Installation

1. Clone repo
```
git clone https://github.com/singleron-RD/CeleScope.git
```

2. Create conda environment and install conda packages. It is recommended to use [mamba](https://github.com/mamba-org/mamba) (which is a faster replacement for Conda)
```
cd CeleScope
conda install mamba
mamba create -n celescope -y --file conda_pkgs.txt
```



3. Install celescope

Make sure you have activated the `celescope` conda environment before running `pip install celescope`. 
```
conda activate celescope
pip install celescope
```

If you want to use `flv_CR` and `flv_trust4`
```
conda activate celescope
pip install celescope>=1.11.0b0
```