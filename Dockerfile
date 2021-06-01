# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

WORKDIR /app
COPY conda_pkgs.txt conda_pkgs.txt
RUN conda install --file conda_pkgs.txt \
 --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing \
 && conda clean --all
COPY requirements.txt requirements.txt
RUN pip install -i https://pypi.tuna.tsinghua.edu.cn/simple -r requirements.txt
COPY . .
RUN python setup.py install

