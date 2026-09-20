# ---------------------------------------------------------------------------
# EDLRM — Dockerfile
#
# Reproducible deployment image for the EDLRM Shiny web server, intended for
# institutional / local hosting as an alternative to the public shinyapps.io
# deployment (see README > Large reference dataset for why this matters).
#
# Build:
#   docker build -t edlrm .
#
# Run (mounting the large reference dataset separately, see README):
#   docker run -p 3838:3838 \
#     -v /path/to/sc1gexpr.h5:/srv/shiny-server/edlrm/sc1gexpr.h5 \
#     edlrm
#
# Then open http://localhost:3838
#
# NOTE: this image has not been exhaustively benchmarked; treat it as a
# reproducibility aid / starting point rather than a production-hardened
# artefact. PRs improving it are welcome.
# ---------------------------------------------------------------------------

FROM rocker/shiny:4.3.2

# System libraries required by hdf5r, Seurat and reticulate/TensorFlow
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libglpk40 \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# R package dependencies
COPY install.R /tmp/install.R
RUN Rscript /tmp/install.R

# Python backend for keras/tensorflow (via reticulate)
RUN Rscript -e 'keras::install_keras(method = "virtualenv")'

# App code and small reference assets (the ~1.2 GB sc1gexpr.h5 backend is
# intentionally excluded — mount it at runtime, see usage note above)
WORKDIR /srv/shiny-server/edlrm
COPY app.R genes_neural.csv my_matrix_filled_merged.csv \
     my_model_MM_prueba_2_all.h5 \
     sc1conf.rds sc1def.rds sc1gene.rds sc1meta.rds ./
COPY www ./www

EXPOSE 3838
CMD ["/usr/bin/shiny-server"]
