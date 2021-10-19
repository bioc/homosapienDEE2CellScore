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
