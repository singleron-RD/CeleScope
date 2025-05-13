FROM mambaorg/micromamba:0.22.0 AS build
ADD . /usr/src/celescope
WORKDIR /usr/src/celescope
USER root
RUN micromamba create --name runtime --always-copy --file conda_pkgs.txt

FROM mambaorg/micromamba:0.22.0 AS runtime
USER root
ADD . /usr/src/celescope
WORKDIR /usr/src/celescope
ENV ENV_NAME=runtime PATH="/opt/conda/envs/runtime/bin:${PATH}"
COPY --from=build --chown=$MAMBA_USER:$MAMBA_USER /opt/conda/envs/runtime /opt/conda/envs/runtime
RUN pip install --verbose --no-cache-dir .
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils less procps && rm -rf /var/lib/apt/lists/*
