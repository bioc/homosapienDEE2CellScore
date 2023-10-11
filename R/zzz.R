#' @importFrom ExperimentHub createHubAccessors
#' @importFrom utils read.csv
.onLoad <- function(libname, pkgname) {
   fl <- system.file("extdata", "metadata.csv", package=pkgname)
   titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
   suppressMessages(createHubAccessors(pkgname, titles))
}


#' Automatically created ergonomic accessor functions
#'
#' Accessor functions for retrieving the data associated with this
#' data package from ExperimentHub. Each of these functions downloads the
#' container file and then returns a path to it. This file can be
#' rehydrated into a SummarizedExperiment by using `readInSEZip`.
#' Usually you would want to actually use `downloadAllTheData` instead of
#' using any of these functions.
#'
#' \itemize{
#'   \item HomosapienDEE2_QC_WARN_Raw    Raw data including data that has quality control warnings
#'   \item HomosapienDEE2_QC_PASS_Raw    Raw data without any quality control warnings
#'   \item HomosapienDEE2_QC_WARN_Rank   Rank normalised data including data that has quality control warnings
#'   \item HomosapienDEE2_QC_PASS_Rank   Rank normalised data without any quality control warnings
#'   \item HomosapienDEE2_QC_WARN_Agg    Aggregated data including data that has quality control warnings
#'   \item HomosapienDEE2_QC_PASS_Agg    Aggregated data without any quality control warnings
#'   \item HomosapienDEE2_QC_WARN_Deseq2 Deseq2 normalised data that has quality control warnings
#'   \item HomosapienDEE2_QC_PASS_Deseq2 Deseq2 normalised data without any quality control warnings
#' }
#'
#' @name HomosapienDEE2_QC_WARN_Raw
#' @aliases HomosapienDEE2_QC_PASS_Raw
#' @aliases HomosapienDEE2_QC_WARN_Rank
#' @aliases HomosapienDEE2_QC_PASS_Rank
#' @aliases HomosapienDEE2_QC_WARN_Agg
#' @aliases HomosapienDEE2_QC_PASS_Agg
#' @aliases HomosapienDEE2_QC_WARN_Deseq2
#' @aliases HomosapienDEE2_QC_PASS_Deseq2
#' @seealso readInSEZip
#' @seealso downloadAllTheData
#' @examples
#' # The ExperimentHub metadata for the Deseq2 normalised data that passes QC is downloadable like so
#' the_metadata <- HomosapienDEE2_QC_PASS_Deseq2(metadata=TRUE)
#' 
#' # Or to download all of the data for the Deseq2 normalised data that passes QC do the following
#' #the_data <- HomosapienDEE2_QC_PASS_Deseq2()
NULL
