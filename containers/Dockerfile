FROM rocker/r-ver:4.3.3

LABEL maintainer="Inkasimo" \
      org.opencontainers.image.title="multiome-link-ranking" \
      org.opencontainers.image.description="Containerized R/Snakemake workflow for multiome peak-gene link ranking benchmarks"

ENV DEBIAN_FRONTEND=noninteractive
ENV RENV_CONFIG_REPOS_OVERRIDE=https://packagemanager.posit.co/cran/2024-01-01
ENV RENV_PATHS_LIBRARY=/opt/renv/library
ENV RENV_PATHS_CACHE=/opt/renv/cache
ENV RENV_CONFIG_CACHE_SYMLINKS=FALSE

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip \
    git curl wget ca-certificates \
    build-essential \
    zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    libhdf5-dev \
    libpng-dev libjpeg-dev libtiff5-dev \
    libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    pkg-config libfontconfig1-dev libcairo2-dev libglib2.0-dev \
    libglpk-dev libgomp1 \
    tabix samtools bedtools \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages('renv', repos=Sys.getenv('RENV_CONFIG_REPOS_OVERRIDE'))"
RUN R -q -e "install.packages('BiocManager', repos=Sys.getenv('RENV_CONFIG_REPOS_OVERRIDE'))"

RUN pip3 install --no-cache-dir \
    snakemake==7.32.4 \
    pulp==2.7.0 \
    "pyyaml>=6.0"

WORKDIR /work

COPY renv.lock /work/renv.lock
COPY renv/activate.R /work/renv/activate.R
COPY .Rprofile /work/.Rprofile

RUN R -q -e '\
renv::settings$bioconductor.version("3.18"); \
renv::restore(project="/work", lockfile="/work/renv.lock", prompt=FALSE, rebuild=TRUE) \
'

CMD ["bash"]
