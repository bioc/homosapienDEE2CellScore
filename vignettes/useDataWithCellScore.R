BiocManager::install()
BiocManager::install("getDEE2")
library(SummarizedExperiment)
library(getDEE2)
library(devtools)
devtools::install_github(repo="flaviusb/CellScore")
library(CellScore)
#load("homosapien_data.RData")
# get data from hub Accessors
library(DESeq2)
sm <- all_built_homosapien_data$qc_pass_deseq2
test1 <- sm[, sm$category == 'test']
standard <- sm[, sm$category == 'standard']
sm1 <- cbind(test1, standard)
start <- c("FIB", "FIB")
test <- c ("iPS-FIB", "iHEP-FIB")
target <- c("ESC", "ESC")
cell.change <- data.frame(start, test, target)
group.OnOff <- OnOff(sm1, cell.change, out.put="marker.list")
