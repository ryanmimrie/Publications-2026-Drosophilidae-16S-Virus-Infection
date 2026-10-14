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

## Repository contents

### `data/`

| File | Description |
| --- | --- |
| `reads/raw/` | Paired-end 16S rRNA sequencing reads (`*_R1.fastq.gz` and `*_R2.fastq.gz`), including negative controls. |
| `data_metadata.csv` | Sample identifiers, host species, experimental conditions and sample-processing metadata. |
| `data_drosophilidae.csv` | Host species names and diet information. |
| `data_viralLoad.csv` | DCV viral-load measurements by host species and experimental block. |
| `data_OTUs.csv` | DADA2 sequence variants, nucleotide sequences and taxonomic assignments. |
| `data_ASVs.csv` | Lookup between ASV identifiers, sequence hashes and labels used in the figures. |
| `data_phyloseq.RData` | Processed phyloseq object used in downstream analyses. |
| `tree_hosts.nwk` | Host phylogeny. |
| `tree_OTUs.nwk` | Bacterial phylogenetic trees. |
| `tree_ASVs.nex` | Tree used to arrange the ASVs in the heatmap. |
| `silva_*` | SILVA v138.1 reference files used for taxonomic assignment; large reference files are provided in compressed parts where applicable. |

### `scripts/`

| Script | Description |
| --- | --- |
| `00_setup.R` | Install R package dependencies. |
| `01_Bioinformatics.R` | Process raw sequencing reads with DADA2, assign taxonomy and prepare the microbiome data. |
| `02_Heatmap.R` | Plot bacterial ASV relative abundances across host species and experimental conditions, alongside phylogenies. |
| `03_AlphaDiversities.R` | Analyse alpha diversity, phylogenetic heritability and repeatability. |
| `04_PCA.R` | Principal component analysis and associated phylogenetic models. |
| `05_NMDS.R` | Non-metric multidimensional scaling and associated phylogenetic models. |
| `06_ASVByAbundances.R` | Model changes in bacterial ASV abundance across conditions. |
| `07_ConditionByAbundance.R` | Model correlations in bacterial abundance between experimental conditions. |
| `08_VirusByAbundance.R` | Model relationships between bacterial abundance and DCV viral load. |

### Other directories

- `models/` — saved R model fits used by the analysis scripts.
- `img/` — repository image assets.
