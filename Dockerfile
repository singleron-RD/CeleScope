FROM continuumio/miniconda3

RUN git clone https://github.com/zhouyiqi91/CeleScope.git \
 && cd CeleScope \
 && conda install --file conda_pkgs.txt --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing \
 && pip install -i https://pypi.tuna.tsinghua.edu.cn/simple celescope \
 && conda clean --all
