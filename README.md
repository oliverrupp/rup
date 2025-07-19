# rup

RNAseq usability pipeline

# Installation

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("getopt", "ggplot2", "reshape2", "pheatmap", "fastqcr", "Rfastp", "Rsubread"))
```

```bash
R -e "fastqcr::fastqc_install()" 
```

# Usage

```bash
Usage: Projects/rup/rup.R [-[-datafolder|d] <character>] [-[-threads|t] <integer>] [-[-help|h]]

Options:
 -d [folder] location of the data folder
 -t [number] number of threads to use
 -h          print this help message
```

# Prepare Files
