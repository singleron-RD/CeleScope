git clone https://github.com/singleron-RD/CeleScope.git

conda create -n celescope
conda activate celescope
conda install --file conda_pkgs.txt --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing

pip install celescope
