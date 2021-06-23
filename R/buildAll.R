# A build is a vector of tasks.
# A task is a function which takes a state object and returns a state object.
# We create tasks which are strung together to build this dataset; some are simple wrappers around external functions.

# Run a build
runBuild <- function(build) {
  accumulators = list(NULL)
  for (x in build) {
    accumulators_new = list()
    for (accumulator in accumulators) {
      for (inner in x) {
        accumulators_new = list(accumulators_new, list(inner(accumulator)))
      }
    }
    accumulators <- accumulators_new
  }
  return(accumulators)
}

# Run a vector of builds - this is just a convenience wrapper around lapply
runBuilds <- function(builds) {
  lapply(builds, runBuild)
}

buildGetData <- function(species, accessions, out_name, metadata=getDEE2Metadata(species)) {
  ret <- function(accumulator) {
    #accumulator[out_name] <- getDEE2::getDEE2(species, accessions)
    #metadata <- getDEE2Metadata(species)
    accumulator[out_name] <- list(do.call(cbind, lapply(accessions, function(y) { getDEE2::getDEE2(species, y, metadata=metadata) })))
    return(accumulator)
  }
  return(ret)
}

buildFilterQC1 <- function() {
  ret <- function(accumulator) {
    accumulator$gene_data <- accumulator$gene_data[accumulator$gene_data$QC_summary == "PASS",]
    return(accumulator)
  }
}

buildFilterQC2 <- function() {
  ret <- function(accumulator) {
    accumulator$gene_data <- accumulator$gene_data[accumulator$gene_data$QC_summary != "FAIL",]
    return(accumulator)
  }
}

cols <- read.csv(system.file("inst", "hsapiens_colData.csv", package="homosapienDEE2CellScore"))
# A vector of the builds that create the `inst` directory are here:
createInst = c(list(c(buildGetData("hsapiens", as.list(cols$SRR_accession[1:10]), "gene_data"))), list(c(buildFilterQC1(), buildFilterQC2())))

