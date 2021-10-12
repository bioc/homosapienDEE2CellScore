
# hsapiens column data we are targetting; notably the accession are in cols$SRR_accession
cols <- read.csv(system.file("inst", "hsapiens_colData.csv", package="homosapienDEE2CellScore"))

#' buildData builds the data included in this package
#'
#' This function generates the data set for this package.
#' All parameters are optional; by default the function will
#' generate a normalised dataset based on downloading
#' the accessions in `inst/hsapiens_colData.csv` for species "hsapiens",
#' and save the dataset to a file called `homosapienDEE2Data.rds` in the current directory.
#'
#' @param species       The species to fetch data for; default is "hsapiens".
#' @param name_prefix   The output file name prefix; default is "homosapienDEE2Data".
#' @param name_suffix   The output file name suffix; default is ".csv"
#' @param build_deseq2  Whether to build the deseq2 normalisation.
#' @param base          The directory to output the file to; default is the current working directory.
#' @param quiet         Whether to suppress notification output where possible; default TRUE.
#' @param metadata      If you have already downloaded metadata for the species, you can pass it in here. Otherwise the metadata will be downloaded.
#' @param counts.cutoff Cutoff value for minimum gene expression; default is 10.
#' @param accessions    Which gene accessions to download data for from DEE2; default is derived from `hsapiens_colData.csv` in this package. For subsets, you can see the internal `cols` objects `SRR_accession` member.
#' @param in_data       If you have already downloaded the accession data from DEE2, you can pass it through here. Otherwise this data will be downloaded.
#' @param dds_design    The design formula used as part of DESeq2 normalisation. Default is `~ 1`. See the documentation for `DESeq2::DESeqDataSetFromMatrix` for more details.
#' @export
#' @import SummarizedExperiment
#' @importFrom getDEE2 getDEE2
#' @importFrom getDEE2 getDEE2Metadata
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom BiocGenerics estimateSizeFactors counts cbind
#' @examples
#' # To build the default, full dataset, and write it out to several csv files:
#' #homosapienDEE2CellScore::buildData()
#'
#' # To build a restricted set of data, with a cached metadata file, only running deseq2 normalisation, to "data_PASS_deseq2.csv" and "data_WARN_deseq2.csv"
#' metadata <- getDEE2Metadata("hsapiens", quiet=TRUE)
#' homosapienDEE2CellScore::buildData(metadata=metadata, accessions=as.list(cols$SRR_accession[1:10]), build_deseq2=TRUE, name_prefix="data")

buildData <- function(species="hsapiens", name_prefix="homosapienDEE2Data", name_suffix=".csv", build_deseq2=TRUE, base=getwd(), quiet=TRUE, metadata=getDEE2Metadata(species, quiet=quiet), counts.cutoff = 10, accessions=as.list(cols$SRR_accession), in_data = do.call(cbind, lapply(accessions, function(y) { getDEE2::getDEE2(species, y, metadata=metadata, quiet=quiet) })), dds_design = ~ 1) {

  out <- list()
  outputs <- list()
  # All of the 'optionally overriden' data possible is calculated in the function's arguments.

  # Take either the clean pass of qc data, or the data which passes and the data which has warnings, but not failures
  qc_pass <- in_data[, startsWith(in_data$QC_summary, "PASS")]
  qc_warn <- in_data[, startsWith(in_data$QC_summary, "PASS") | startsWith(in_data$QC_summary, "WARN")]

  # Now we filter based on gene activity
  qc_pass_filtered <- qc_pass[rowSums(assay(qc_pass)) > counts.cutoff,]
  qc_warn_filtered <- qc_warn[rowSums(assay(qc_warn)) > counts.cutoff,]

  if(build_deseq2) {
    # Now normalisation
    dds_qc_pass_filtered <- BiocGenerics::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(
      countData = SummarizedExperiment::assay(qc_pass_filtered, "counts"),
      colData = SummarizedExperiment::colData(qc_pass_filtered),
      rowData = SummarizedExperiment::rowData(qc_pass_filtered),
      design = dds_design))
    logcounts_qc_pass_filtered <- log2(counts(dds_qc_pass_filtered, normalize=TRUE) + 1)
    dds_qc_warn_filtered <- BiocGenerics::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(
      countData = SummarizedExperiment::assay(qc_warn_filtered, "counts"),
      colData = SummarizedExperiment::colData(qc_warn_filtered),
      rowData = SummarizedExperiment::rowData(qc_warn_filtered),
      design = dds_design))
    logcounts_qc_warn_filtered <- log2(counts(dds_qc_warn_filtered, normalize=TRUE) + 1)
    out <- c(out, list(qc_pass_deseq2=logcounts_qc_pass_filtered, qc_warn_deseq2=logcounts_qc_warn_filtered))
    outputs <- c(outputs, list(qc_pass_deseq2=paste(name_prefix, "_PASS_deseq2", name_suffix, sep=""), qc_warn_deseq2=paste(name_prefix, "_WARN_deseq2", name_suffix, sep="")))
  }

  ## And beginning of pca
  ##pca_qc_pass_filtered <- prcomp(t(logcounts_qc_pass_filtered))
  ##pca_qc_warn_filtered <- prcomp(t(logcounts_qc_warn_filtered))
  
  # Write results out, gathered into a named list
  #out <- list(qc_pass_deseq2=logcounts_qc_pass_filtered, qc_pass_deseq2_pca=pca_qc_pass_filtered, qc_warn_deseq2=logcounts_qc_warn_filtered, qc_warn_deseq2_pca=pca_qc_warn_filtered)
  #saveRDS(out, file=paste(base, name, sep="/"))
  writeOutput(out, outputs=outputs)
}

# Takes a named-list of processed data plus a bunch of filenames for each named thing
writeOutput <- function(the_data, outputs=list(dds_qc_pass_filtered="dds_qc_pass_filtered.csv", dds_qc_warn_filter="dds_qc_warn_filter.csv")) {
  lapply(names(outputs), function(n) {
    write.csv(the_data[[n]], file=outputs[[n]], row.names=TRUE)
  })
}
