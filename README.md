# homosapienDEE2CellScore
Example packaged cell data

## Installation

This package is not yet uploaded to Bioconductor; in order to install it from this git repo, first clone the repo:

```sh
git clone https://github.com/flaviusb/homosapienDEE2CellScore.git
```

Then, run R inside the cloned directory. In order to load the package locally, you can run:

```R
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("getDEE2")
devtools::load_all()
```

## Generating the packaged data from HEAD

To generate the data from the version of this package on github, you can run the following in R:

```R
BiocManager::install()
BiocManager::install("getDEE2")
library(Biobase)
library(SummarizedExperiment)
library(getDEE2)
library(devtools)
devtools::install_github(repo="flaviusb/CellScore")
devtools::install_github(repo="flaviusb/homosapienDEE2CellScore")
library(homosapienDEE2CellScore)
buildData(build_tsne=FALSE)
```
