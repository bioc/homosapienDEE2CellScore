# A build is a list of bundles, a bundle is a list of tasks.
# We accumulate the cross product of the tasks in each bundle
# A task is a function which takes a state object and returns a state object.
# We create tasks which are strung together to build this dataset; some are simple wrappers around external functions.

# Run a build
runBuild <- function(build) {
  accumulators = list(list(tags=list("initial", "padding"), padding=list("padding2", "padding3")))
  for (x in build) {
    accumulators_new = list()
    for (accumulator in accumulators) {
      print("Accumulator at beginning of for loop: ")
      print(accumulator)
      for (inner in x) {
        accumulators_new = c(accumulators_new, (inner(accumulator)))
      }
    }
    accumulators <- accumulators_new
  }
  return(accumulators)
}

# Run a list of builds - this is just a convenience wrapper around lapply
runBuilds <- function(builds) {
  lapply(builds, runBuild)
}

buildGetData <- function(species, accessions, out_name, metadata=getDEE2Metadata(species, quiet=TRUE)) {
  ret <- function(accumulator) {
    print("accumulator at start of getData: ")
    print(accumulator)
    accumulator$tags <- c(accumulator$tags, "foo")
    accumulator[[out_name]] <- list(do.call(cbind, lapply(accessions, function(y) { getDEE2::getDEE2(species, y, metadata=metadata, quiet=TRUE) })))
    print("Accumulator at end of getData: ")
    print(accumulator)
    return(accumulator)
  }
  return(ret)
}

buildFilter <- function(filt, on, tag) {
  ret <- function(accumulator) {
    print("accumulator as of build filter")
    print(accumulator)
    r <- accumulator[on]
    print("r: ")
    print(r)
    m <- list(r[,filt(r)])

    accumulator[on] <- m
    accumulator$tag <- c(accumulator$tag, tag)
    return(accumulator)
  }
  return(ret)
}

filtQC1 <- buildFilter(function(it) { startsWith(it$QC_summary, "PASS") }, "gene_data", "filter: pass")
filtQC2 <- buildFilter(function(it) { startsWith(it$QC_summary, "PASS") | startsWith(it$QC_summary, "WARN") }, "gene_data", "filter: pass and warn")
filtNoQC <- buildFilter(function(it) { it$QC_summary != "TEST" }, "gene_data", "filter: no filter")

# add c(printAccumulator) to a point in createInst to see what the accumulators are there
mkPrintAccumulator <- function(message) {
  foo <- function(accumulator) {
    printAccumulator(accumulator, message=message)
  }
  foo
}
printAccumulator <- function(accumulator, message="") {
  if(nchar(message) > 0) {
    print(message)
  }
  print(accumulator)
  accumulator
}

cols <- read.csv(system.file("inst", "hsapiens_colData.csv", package="homosapienDEE2CellScore"))
# A list of the builds that create the `inst` directory are here:
createInst = list(
  c(mkPrintAccumulator(message="Initial accumulator:")),
  c(buildGetData("hsapiens", as.list(cols$SRR_accession[285:295]), "gene_data")),
  c(mkPrintAccumulator(message="Accumulator after buildGetData:")),
  (c(filtQC1, filtQC2, filtNoQC)),
  c(mkPrintAccumulator(message="Final accumulator:")))

doIt <- function() {
  species <- "hsapiens"
  accessions <- as.list(cols$SRR_accession[285:295])
  metadata <- getDEE2Metadata(species, quiet=TRUE)
  in_data <- do.call(cbind, lapply(accessions, function(y) { getDEE2::getDEE2(species, y, metadata=metadata, quiet=TRUE) }))
  #print(in_data)
  qc_pass <- in_data[, startsWith(in_data$QC_summary, "PASS")]
  qc_warn <- in_data[, startsWith(in_data$QC_summary, "PASS") | startsWith(in_data$QC_summary, "WARN")]
  list(qc_pass=qc_pass, qc_warn=qc_warn)
}
