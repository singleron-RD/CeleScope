# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

WORKDIR /app
COPY conda_pkgs.txt conda_pkgs.txt
RUN mkdir -p /opt/conda/pkgs/cache && conda clean --packages && conda clean --all \
&& conda install -c conda-forge mamba &&  mamba install -y --file conda_pkgs.txt && pip install celescope

