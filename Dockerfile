# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

WORKDIR /app
COPY conda_pkgs.txt conda_pkgs.txt
RUN conda install -c conda-forge mamba
RUN mamba install -y --file conda_pkgs.txt
RUN pip install -i https://pypi.mirrors.ustc.edu.cn/simple/ celescope

