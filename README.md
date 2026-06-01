<img src="img/logo.png" width="200" alt="Logo" style="display: block; margin: auto;">

# (2026) Examining interactions between the microbiome and viral infection among Drosophila species

This repository contains data, code, and models used in Imrie et al., (2026) "Examining interactions between the microbiome and viral infection among Drosophila species."

## Quick Start

### Platform Support
![Windows](https://img.shields.io/badge/Windows-blue?logo=microsoftwindows)
![macOS](https://img.shields.io/badge/macOS-black?logo=apple)
![Linux](https://img.shields.io/badge/Linux-grey?logo=linux)

### Dependencies
![R Version](https://img.shields.io/badge/R-4.4.2-blue)
![devtools](https://img.shields.io/badge/devtools-2.4.5-ff69b4)
![tidyverse](https://img.shields.io/badge/tidyverse-2.0.0-blue)

### System Requirements
![RAM](https://img.shields.io/badge/minimum%20RAM-8GB-important)

The included script `scripts/00_setup.R` can be used to install all package dependencies at once.

Scripts in this repository use the `here` library to dynamically set paths. For this to work correctly, Rstudio must be opened by double-clicking on one of the files in `scripts`. Path errors will appear if Rstudio was first opened using a shortcut or a script from a different location.

## Contents

| Directory               | Description                                                                             |
|-------------------------|-----------------------------------------------------------------------------------------|
| `data/`                 | Contains all data files used in this study                                              |
| └─ `reads`              | Contains raw and filtered read files in fastq.gz format                                 |
| `img/`                  | Image files used in this repository                                                     |
| `models/`               | Contains all MCMCglmm files analysed in this study                                      |
| `scripts/`              | Contains all scripts used in this study including bioinformatics, modelling, and plots  |
| `.here`                 | Empty text file used for dynamic pathing by the here() library                          |
| `LICENCE.txt`           | CC-BY-4.0 license designation                                                           |
| `README.md`             | Markdown file for repository documentation                                              |

