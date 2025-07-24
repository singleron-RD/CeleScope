# Build stage
FROM mambaorg/micromamba:0.22.0 AS build
ADD . /usr/src/celescope
WORKDIR /usr/src/celescope
USER root

# conda pkgs
RUN micromamba create --name runtime --always-copy --file conda_pkgs.txt \
    && micromamba clean --all --yes

# Runtime stage
FROM mambaorg/micromamba:0.22.0 AS runtime
ADD . /usr/src/celescope
WORKDIR /usr/src/celescope
USER root
ENV ENV_NAME=runtime PATH="/opt/conda/envs/runtime/bin:${PATH}"

# Copy the environment from the build stage
COPY --from=build --chown=$MAMBA_USER:$MAMBA_USER /opt/conda/envs/runtime /opt/conda/envs/runtime

# Install Python dependencies and necessary system packages in one layer
RUN pip install --verbose --no-cache-dir . \
    && apt-get update \
    && apt-get install -y --no-install-recommends apt-utils less procps \
    && rm -rf /var/lib/apt/lists/*

# Optionally, switch back to a non-root user if desired
USER $MAMBA_USER

