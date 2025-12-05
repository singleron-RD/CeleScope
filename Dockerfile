# =========================
# Build stage
# =========================
FROM mambaorg/micromamba:0.22.0 AS build

# Set working directory
WORKDIR /usr/src/celescope

# Copy source code
COPY . .

# Switch to root to install conda packages
USER root

# Create conda environment and clean cache
RUN micromamba create --name runtime --always-copy --file conda_pkgs.txt \
    && micromamba clean --all --yes

# Install Python dependencies in the conda environment
RUN micromamba run -n runtime pip install --verbose --no-cache-dir .

# =========================
# Runtime stage
# =========================
FROM mambaorg/micromamba:0.22.0 AS runtime

WORKDIR /usr/src/celescope

# Copy only the installed environment from build stage
COPY --from=build /opt/conda/envs/runtime /opt/conda/envs/runtime

# Set environment variables
ENV ENV_NAME=runtime
ENV PATH="/opt/conda/envs/runtime/bin:${PATH}"

# Install minimal system dependencies
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends apt-utils less procps \
    && rm -rf /var/lib/apt/lists/*

# Switch back to non-root user
USER mambauser
