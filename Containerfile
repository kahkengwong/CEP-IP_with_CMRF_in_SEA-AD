# Containerfile for the SEA-AD CEP-IP/CMRF analysis
#
# R 4.6.1 with the packages in renv.lock, and Python 3.13.14 with the
# packages in requirements-python.txt, as used for the paper.
#
# Build (from the folder that holds this file, scripts/, verification/, seaad-setup.sh,
# seaad-run.sh, renv.lock and requirements-python.txt):
#     podman build --no-cache -t seaad-cepip -f Containerfile .
#
# Run the whole analysis (the SEA-AD data folder is mounted at /data):
#     podman run --rm --shm-size=16g -v /path/to/SEA-AD:/data:Z seaad-cepip seaad-run all
#
# See README_container.md for the inputs, the run modes and the run order.

FROM docker.io/rocker/r-ver:4.6.1

LABEL org.opencontainers.image.title="SEA-AD CEP-IP/CMRF analysis" \
      org.opencontainers.image.source="https://github.com/kahkengwong"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# ---- System libraries ------------------------------------------------------
# Build tools, and the libraries the R packages compile against:
#   sf, s2, units          GDAL, GEOS, PROJ, udunits, sqlite, abseil (monocle3)
#   BPCells, rhdf5         HDF5
#   ragg, Cairo, systemfonts, textshaping   graphics and font libraries
#   gsl, igraph, nloptr    GSL, GLPK, NLopt
#   htmlwidgets            pandoc, for the self-contained 3-D HTML figures
# fonts-liberation provides Liberation Sans and Liberation Sans Narrow, which
# fontconfig uses in place of Arial and Arial Narrow (same letter widths).
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl git bzip2 pkg-config cmake make \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libgdal-dev libgeos-dev libproj-dev libudunits2-dev libsqlite3-dev \
        libabsl-dev libhdf5-dev \
        libcairo2-dev libxt-dev libfreetype6-dev libpng-dev libtiff-dev \
        libjpeg-dev libwebp-dev libharfbuzz-dev libfribidi-dev \
        libfontconfig1-dev \
        libgsl-dev libglpk-dev libnlopt-dev libgmp-dev \
        zlib1g-dev libbz2-dev liblzma-dev libzstd-dev libicu-dev \
        pandoc fontconfig fonts-liberation \
    && rm -rf /var/lib/apt/lists/* \
    && fc-cache -f

# ---- R packages ------------------------------------------------------------
# renv installs the exact versions in renv.lock into the site library.
# CRAN packages come from Posit Package Manager as Linux binaries (faster);
# Bioconductor 3.23 packages and the three GitHub packages (monocle3, presto,
# BPCells) are built from source at the commits recorded in the lock file.
# The lock file records pals and svglite as absent; the scripts fall back to
# base R for both, as in the runs for the paper.
WORKDIR /opt/seaad
COPY renv.lock /opt/seaad/renv.lock
# renv installs one package at a time (RENV_CONFIG_INSTALL_JOBS=1). The parallel installation of
# renv >= 1.2 can load-test a package while one of its dependencies is still being installed.
RUN . /etc/os-release \
    && export RENV_CONFIG_REPOS_OVERRIDE="https://packagemanager.posit.co/cran/__linux__/${VERSION_CODENAME}/latest" \
    && export RENV_CONFIG_INSTALL_JOBS=1 \
    && Rscript -e 'install.packages("renv", repos = Sys.getenv("RENV_CONFIG_REPOS_OVERRIDE"))' \
    && Rscript -e 'options(Ncpus = max(1L, parallel::detectCores() - 1L)); \
                   renv::restore(project = "/opt/seaad", lockfile = "/opt/seaad/renv.lock", \
                                 library = .libPaths()[1], prompt = FALSE)' \
    && rm -rf /opt/seaad/renv /tmp/Rtmp* /root/.cache/R

# ---- Python ----------------------------------------------------------------
# Python 3.13.14 from conda-forge (as on the Fedora workstation), then the
# pinned packages with pip.
COPY requirements-python.txt /opt/seaad/requirements-python.txt
RUN curl -fsSL -o /tmp/miniforge.sh \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh \
    && /opt/conda/bin/conda create -y -n seaad -c conda-forge python=3.13.14 \
    && /opt/conda/envs/seaad/bin/python -m pip install --no-cache-dir \
           -r /opt/seaad/requirements-python.txt \
    && /opt/conda/bin/conda clean -afy \
    && ln -s /opt/conda/envs/seaad/bin/python /usr/local/bin/python \
    && ln -s /opt/conda/envs/seaad/bin/python /usr/local/bin/python3

# ---- Scripts ---------------------------------------------------------------
# The scripts find the project folder through SEAAD_ROOT (the mounted /data).
ENV SEAAD_ROOT=/data
COPY scripts/R/ /opt/seaad/R/
COPY scripts/Python/ /opt/seaad/Python/
COPY scripts/SuppTables/ /opt/seaad/SuppTables/
COPY verification/compare_rerun.py verification/compare_python.py verification/compare_gamdata.py \
     /opt/seaad/verification/
COPY seaad-setup.sh /usr/local/bin/seaad-setup
COPY seaad-run.sh /usr/local/bin/seaad-run
RUN chmod +x /usr/local/bin/seaad-setup /usr/local/bin/seaad-run

# Versions, for the record (written at build time)
RUN { echo "R:"; R --version | head -1; \
      echo "Python:"; python --version; \
      Rscript -e 'ip <- installed.packages()[, c("Package", "Version")]; write.table(ip, sep = "\t", quote = FALSE, row.names = FALSE)'; \
      python -m pip freeze; } > /opt/seaad/versions.txt 2>&1

WORKDIR /data
CMD ["bash"]
