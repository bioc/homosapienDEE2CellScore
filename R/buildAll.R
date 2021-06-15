# A build is a vector of tasks.
# A task is a function which takes a state object and returns a state object.
# We create tasks which are strung together to build this dataset; some are simple wrappers around external functions.
# A vector of the builds that create the `inst` directory are here:
createInst = c()

# Run a build
runBuild <- function(build) {
  accumulator = NULL
  for (x in build) {
    accumulator = call(x, accumulator)
  }
}

# Run a vector of builds - this is just a convenience wrapper around lapply
runBuilds <- function(builds) {
  lapply(builds, runBuild)
}
