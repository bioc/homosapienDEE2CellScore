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
        accumulators_new = c(accumulators_new, (inner(accumulator)))
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
    accumulator[out_name] <- list(do.call(cbind, lapply(accessions, function(y) { getDEE2::getDEE2(species, y, metadata=metadata) })))
    return(accumulator)
  }
  return(ret)
}

buildFilter <- function(filt, on) {
  ret <- function(accumulator) {
    m <- accumulator[,filt(accumulator)]

    accumulator <- m
    return(accumulator)
  }
  return(ret)
}

filtQC1 <- buildFilter(function(it) { startsWith(it$QC_summary, "PASS") }, "gene_data")
filtQC2 <- buildFilter(function(it) { startsWith(it$QC_summary, "PASS") | startsWith(it$QC_summary, "WARN") }, "gene_data")
filtNoQC <- buildFilter(function(it) { it$QC_summary != "TEST" }, "gene_data")


cols <- read.csv(system.file("inst", "hsapiens_colData.csv", package="homosapienDEE2CellScore"))
# A vector of the builds that create the `inst` directory are here:
createInst = list(c(buildGetData("hsapiens", as.list(cols$SRR_accession[285:295]), "gene_data")), (c(filtQC1, filtQC2, filtNoQC)))

