#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# EDLRM — dependency installation script
#
# Installs every CRAN package app.R needs. Run once before launching the app
# locally:
#
#   Rscript install.R
#
# keras/tensorflow additionally need a Python backend the first time only:
#
#   Rscript -e 'keras::install_keras()'
# ---------------------------------------------------------------------------

cran_pkgs <- c(
  "shiny", "shinydashboard", "DT", "magrittr", "Matrix", "shinyhelper",
  "data.table", "tensorflow", "reticulate", "hdf5r", "ggdendro",
  "gridExtra", "shinythemes", "ggplot2", "ggrepel", "bslib", "shinyjs",
  "rmarkdown", "shinyalert", "stringr", "isoband", "keras", "Seurat"
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_pkgs, installed)

if (length(to_install) > 0) {
  message("Installing ", length(to_install), " missing package(s): ",
          paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required CRAN packages are already installed.")
}

message("\nNext step (first run only): Rscript -e 'keras::install_keras()'")
