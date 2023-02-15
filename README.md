# homosapienDEE2CellScore
Example packaged cell data

## Requirements

This package requires Bioconductor version 3.17 or later, as that version is when the data of this data package was added to ExperimentHub. That in turn requires R version 4.3.0.

## Installation

This package is not yet uploaded to Bioconductor; in order to install it from this git repo, first clone the repo:

```sh
git clone https://github.com/flaviusb/homosapienDEE2CellScore.git
```

Then, run R inside the cloned directory. In order to load the package locally, you need to set things up like so:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("devtools")
BiocManager::install()
BiocManager::install("getDEE2")
BiocManager::install("DESeq2")
library(DESeq2)
library(Biobase)
library(SummarizedExperiment)
library(getDEE2)
library(devtools)
devtools::install_github(repo="flaviusb/CellScore")
library(CellScore)
```

And then you can run the following to load the package locally:
```R
devtools::load_all()
```

## Generating the packaged data from HEAD

To generate the data from the version of this package on github, you can run the following in R:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
install.packages("devtools")
BiocManager::install()
BiocManager::install("getDEE2")
BiocManager::install("DESeq2")
library(DESeq2)
library(Biobase)
library(SummarizedExperiment)
library(getDEE2)
library(devtools)
devtools::install_github(repo="flaviusb/CellScore")
devtools::install_github(repo="flaviusb/homosapienDEE2CellScore")
library(CellScore)
library(homosapienDEE2CellScore)
buildData(build_tsne=FALSE)
```

## Just getting all the data

There is a helper to download the pregenerated data and format it into a tagged list of SummarizedExperiments. If you already have this library installed, you can use it like so:

```R

the_pregenerated_data <- downloadAllTheData()
```
