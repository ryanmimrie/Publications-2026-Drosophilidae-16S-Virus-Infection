# ==============================================================================
# ===== (2026) Drosophilidae 16S Virus Infection: Dependencies Setup ===========
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script installs all dependencies used in this study.

# Note: This script will install the most recent versions of each dependency.
#       The exact versions of each tool used in this study are provided below.
#       It may be necessary to install exact versions if significant changes
#       have been made to a library since this study was published.

# ------------------------------------------------------------------------------
# ----- 1. Install Packages ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Specify Packages --------------------------------------------------

packages <- data.frame(
  package = c("tidyverse", "phyloseq", "MCMCglmm", "dada2",  "Biostrings", "vegan", "microbiome", "david-barnett/microViz", "zCompositions", "compositions", "ape",   "ggtree", "Cairo", "here",  "progress", "openssl", "conflicted", "decontam", "patchwork", "RColorBrewer", "plotrix", "scales"),
  version = c("2.0.0",     "1.50.0",   "2.36",     "1.34.0", "2.74.1",     "2.7.1", "1.28.0",     "0.12.7",                 "1.5.0.5",       "2.0.8",        "5.8-1", "3.14.0", "1.6.2", "1.0.1", "1.2.3",    "2.3.3",   "1.2.0",      "1.26.0",   "1.3.1",     "1.1.3",        "3.8.4",   "1.4.0"),
  source = c("cran",       "bioc",     "cran",     "bioc",   "bioc",       "cran",  "bioc",       "github",                 "cran",          "cran",         "cran",  "bioc",   "cran",  "cran",  "cran",     "cran",    "bioc",       "bioc",     "cran",      "cran",         "cran",    "cran"))

# ----- 1.2. Install Functions -------------------------------------------------

install_cran <- function(package) {
    install.packages(package, dependencies = TRUE)
}

install_bioc <- function(package) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
  BiocManager::install(package, ask = FALSE, update = FALSE)
}

install_github <- function(repo) {
  if (!requireNamespace("remotes", quietly = TRUE)) {install.packages("remotes")}
  remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
}


# ----- 1.3. Install Packages --------------------------------------------------

for (i in c(1:nrow(packages))) {
  
  package <- packages$package[i]
  source <- packages$source[i]
  
  if (requireNamespace(package, quietly = TRUE)) {
    print(sprintf("%s already installed, skipping.", package))
    next
  }
  
  print(sprintf("Installing %s", package))
  
  if (source == "cran") {
    install_cran(package)
  } else if (source == "bioc") {
    install_bioc(package)
  } else {
    install_github(package)
  }
  
}
