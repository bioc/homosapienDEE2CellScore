# A build is a vector of tasks.
# A task is a function which takes a state object and returns a state object.
# We create tasks which are strung together to build this dataset; some are simple wrappers around external functions.
# A vector of the builds that create the `inst` directory are here:
createInst = c()

# Run a build
function runBuild(build) {
  accumulator = NULL
  
}

# Run a vector of builds - this is just a convenience wrapper around lapply
function runBuilds(builds) {
  lapply(builds, runBuild)
}
