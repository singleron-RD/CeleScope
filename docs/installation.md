# Hardware/Software Requirements

- minimum 32GB RAM(to run STAR aligner)
- conda(or micromamba)
- git

# Installation

1. Clone repo
```
git clone https://github.com/singleron-RD/CeleScope.git
```

2. Create conda environment and install conda packages. 
If you have conda installed, it is recommended to use [mamba](https://github.com/mamba-org/mamba) (which is a faster replacement for Conda):
```
conda install mamba
cd CeleScope
mamba create -n celescope -y --file conda_pkgs.txt
```
Or you can use [micromamba](https://mamba.readthedocs.io/en/latest/installation.html)
```
curl micro.mamba.pm/install.sh | bash
cd CeleScope
micromamba create -n celescope -y --file conda_pkgs.txt
```


3. Install celescope

Make sure you have activated the `celescope` conda environment before running `pip install celescope`. 
```
conda activate celescope
pip install celescope
```