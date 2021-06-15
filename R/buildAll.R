# A build is a vector of tasks.
# A task is a function which takes a state object and returns a state object.
# We create tasks which are strung together to build this dataset; some are simple wrappers around external functions.

# Run a build
runBuild <- function(build) {
  accumulator = NULL
  for (x in build) {
    accumulator = x(accumulator)
  }
}

# Run a vector of builds - this is just a convenience wrapper around lapply
runBuilds <- function(builds) {
  lapply(builds, runBuild)
}

buildGetData <- function(species, accessions, out_name) {
  ret <- function(accumulator) {
    accumulator[out_name] <- getDEE2::getDEE2(species, accessions)
  }
  return(ret)
}

# A vector of the builds that create the `inst` directory are here:
createInst = c(buildGetData("hsapiens", c(), "gene_data"))

