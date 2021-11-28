
# hsapiens column data we are targetting; notably the accession are in cols$SRR_accession
cols <- read.csv(system.file("inst", "hsapiens_colData.csv", package="homosapienDEE2CellScore"))

#' buildRaw gets the raw data in SummarizedExperiment format
#'
#' This function gets the raw Data from dee2 and packages it in a SummarizedExperiment
#'
#' @param species          The species to fetch data for; default is "hsapiens".
#' @param quiet            Whether to suppress notification output where possible; default TRUE.
#' @param metadata         If you have already downloaded metadata for the species, you can pass it in here. Otherwise the metadata will be downloaded.
#' @param accessions       Which sample ids to download from DEE2 (we refer to these as accessions); default is derived from `hsapiens_colData.csv` in this package. For subsets, you can see the internal `cols` objects `SRR_accession` member.
#' @export
#' @import SummarizedExperiment
#' @importFrom getDEE2 getDEE2
#' @importFrom getDEE2 getDEE2Metadata

buildRaw <- function(species="hsapiens", accessions=as.list(cols$SRR_accession), quiet=TRUE, metadata=getDEE2Metadata(species, quiet=quiet)) {
  return(do.call(cbind, lapply(accessions, function(y) { getDEE2::getDEE2(species, y, metadata=metadata, quiet=quiet) })));
}

#' buildData builds the data included in this package
#'
#' This function generates the data set for this package.
#' All parameters are optional; by default the function will
#' generate a normalised dataset based on downloading
#' the accessions in `inst/hsapiens_colData.csv` for species "hsapiens",
#' and save the dataset to a file called `homosapienDEE2Data.rds` in the current directory.
#'
#' @param species          The species to fetch data for; default is "hsapiens".
#' @param name_prefix      The output file name prefix; default is "homosapienDEE2Data".
#' @param name_suffix      The output file name suffix; default is ".csv"
#' @param generate_qc_pass Generate output from the input data that passed quality control
#' @param generate_qc_warn Generate output from the the conbination of input data that passed quality control and input data that had warnings in quality control
#' @param build_deseq2     Whether to build the deseq2 normalisation.
#' @param base             The directory to output the file to; default is the current working directory.
#' @param quiet            Whether to suppress notification output where possible; default TRUE.
#' @param metadata         If you have already downloaded metadata for the species, you can pass it in here. Otherwise the metadata will be downloaded.
#' @param counts.cutoff    Cutoff value for minimum gene expression; default is 10.
#' @param accessions       Which sample ids to download from DEE2 (we refer to these as accessions); default is derived from `hsapiens_colData.csv` in this package. For subsets, you can see the internal `cols` objects `SRR_accession` member.
#' @param in_data          If you have already downloaded the accession data from DEE2, you can pass it through here. Otherwise this data will be downloaded.
#' @param dds_design       The design formula used as part of DESeq2 normalisation. Default is `~ 1`. See the documentation for `DESeq2::DESeqDataSetFromMatrix` for more details.
#' @param write_files      Write out normalised data to files. If this is false, the function will not write out the normalised data, but will only return it.
#' @export
#' @import SummarizedExperiment
#' @importFrom getDEE2 getDEE2
#' @importFrom getDEE2 getDEE2Metadata
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData as.data.frame
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom Rtsne Rtsne
#' @importFrom BiocGenerics estimateSizeFactors counts cbind
#' @examples
#' # To build the default, full dataset, and write it out to several csv files:
#' #homosapienDEE2CellScore::buildData()
#'
#' # To build a restricted set of data, with a cached metadata file, only running deseq2 normalisation, to "data_PASS_deseq2.csv" and "data_WARN_deseq2.csv"
#' metadata <- getDEE2Metadata("hsapiens", quiet=TRUE)
#' homosapienDEE2CellScore::buildData(metadata=metadata, accessions=as.list(cols$SRR_accession[1:10]), build_deseq2=TRUE, name_prefix="data")
#'
#' # Process a subset of the data, but do not write it out into files
#' processed_data <- homosapienDEE2CellScore::buildData(metadata=metadata, accessions=as.list(cols$SRR_accession[1:10]), write_files=FALSE)
#'
#' # Get PCA form of the deseq2 normalised data that passed quality control
#' pca_form <- prcomp(t(processed_data$qc_pass_deseq2))

buildData <- function(species="hsapiens", name_prefix="homosapienDEE2Data", name_suffix=".csv", build_deseq2=TRUE, build_tsne=TRUE, generate_qc_pass = TRUE, generate_qc_warn = TRUE, base=getwd(), quiet=TRUE, metadata=if((!build_deseq2) && (!qc_pass || !qc_warn)) { return(list()); } else { getDEE2Metadata(species, quiet=quiet) }, counts.cutoff = 10, accessions=as.list(cols$SRR_accession), in_data = if((!build_deseq2) && (!qc_pass || !qc_warn)) { return(list()); } else { buildRaw(species=species, accessions=accessions, quiet=quiet, metadata=metadata) }, dds_design = ~ 1, write_files = TRUE) {

  out <- list()
  outputs <- list()
  # Check whether we are not going to do something
  if((!build_deseq2) && (!qc_pass || !qc_warn)) {
    return(out);
  }
  # All of the 'optionally overriden' data possible is calculated in the function's arguments.

  # Take either the clean pass of qc data, or the data which passes and the data which has warnings, but not failures
  qc_pass <- in_data[, startsWith(in_data$QC_summary, "PASS")]
  qc_warn <- in_data[, startsWith(in_data$QC_summary, "PASS") | startsWith(in_data$QC_summary, "WARN")]

  # Now we filter based on gene activity
  qc_pass_filtered <- qc_pass[rowSums(assay(qc_pass)) > counts.cutoff,]
  qc_warn_filtered <- qc_warn[rowSums(assay(qc_warn)) > counts.cutoff,]

  if(build_deseq2) {
    # Now normalisation
    if(generate_qc_pass) {
      dds_qc_pass_filtered <- BiocGenerics::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(qc_pass_filtered, "counts"),
        colData = SummarizedExperiment::colData(qc_pass_filtered),
        rowData = SummarizedExperiment::rowData(qc_pass_filtered),
        design = dds_design))
      logcounts_qc_pass_filtered <- log2(counts(dds_qc_pass_filtered, normalize=TRUE) + 1)
      out <- c(out, list(qc_pass_deseq2=logcounts_qc_pass_filtered))
      outputs <- c(outputs, list(qc_pass_deseq2=paste(name_prefix, "_PASS_deseq2", name_suffix, sep="")))
    }
    if(generate_qc_warn) {
      dds_qc_warn_filtered <- BiocGenerics::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(qc_warn_filtered, "counts"),
        colData = SummarizedExperiment::colData(qc_warn_filtered),
        rowData = SummarizedExperiment::rowData(qc_warn_filtered),
        design = dds_design))
      logcounts_qc_warn_filtered <- log2(counts(dds_qc_warn_filtered, normalize=TRUE) + 1)
      out <- c(out, list(qc_warn_deseq2=logcounts_qc_warn_filtered))
      outputs <- c(outputs, list(qc_warn_deseq2=paste(name_prefix, "_WARN_deseq2", name_suffix, sep="")))
    }
  }
  if(build_tsne) {
    if(generate_qc_pass) {
      tsne_qc_pass_filtered <- Rtsne(unique(t(as.matrix(SummarizedExperiment::as.data.frame(SummarizedExperiment::assay(qc_pass_filtered, "counts"))))), check_duplicates=FALSE)$Y
      out <- c(out, list(qc_pass_tsne=tsne_qc_pass_filtered))
      outputs <- c(outputs, list(qc_pass_tsne=paste(name_prefix, "_PASS_tsne", name_suffix, sep="")))
    }
    if(generate_qc_warn) {
      tsne_qc_warn_filtered <- Rtsne(unique(t(as.matrix(SummarizedExperiment::as.data.frame(SummarizedExperiment::assay(qc_warn_filtered, "counts"))))), check_duplicates=FALSE)$Y
      out <- c(out, list(qc_warn_tsne=tsne_qc_warn_filtered))
      outputs <- c(outputs, list(qc_warn_tsne=paste(name_prefix, "_WARN_tsne", name_suffix, sep="")))
    }
  }

  if (write_files) {
    writeOutput(out, outputs=outputs)
  }
  out
}

# Takes a named-list of processed data plus a bunch of filenames for each named thing
writeOutput <- function(the_data, outputs=list(dds_qc_pass_filtered="dds_qc_pass_filtered.csv", dds_qc_warn_filter="dds_qc_warn_filter.csv")) {
  lapply(names(outputs), function(n) {
    write.csv(the_data[[n]], file=outputs[[n]], row.names=TRUE)
  })
}

#' srx_agg_se is a version of srx_agg that works on SummarizedExperiments
#'
#' This function aggregates runs that represent the same SRA experiment, and reorganises
#' the coldata in the SummarizedExperiment to to be grouped by SRA experiment in order to
#' preserve necessary SummarizedExperiment internal invariants.
#'
#' @param x          A SummarizedExperiment.
#' @param counts     What kind of count; "GeneCounts" for STAR based gene counts, "TxCounts" for kallisto transcript level counts or "Tx2Gene" for transcript counts aggregated to gene level. Default is "GeneCounts"
#' @export
#' @import SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData as.data.frame

srx_agg_se <- function(x,counts="GeneCounts") {
    mds<-colData(x)
    base<-SummarizedExperiment::assay(x, "counts")
    n=nrow(base)
    srx_es<-unique(mds[["SRX_accession"]])
    SRX_assay <- vapply(X=srx_es, function(srx) {
        srrs<-rownames(mds)[which(mds[["SRX_accession"]] %in% srx)]
        if (length(srrs)>1) {
            rowSums(SummarizedExperiment::assay(x, "counts")[,srrs])
        } else {
            SummarizedExperiment::assay(x, "counts")[,srrs]
        }
    } , numeric(n))
    m<-length(srx_es)
    SRX_cols <- sapply(X=srx_es, function(srx) {
        srrs<-rownames(mds)[which(mds[["SRX_accession"]] %in% srx)]
        list(SummarizedExperiment::assay(x, "counts")[,srrs])
    })
    SRX_coldata <- DataFrame(matrix(SRX_cols, dimnames=list(srx_es)))
    return(SummarizedExperiment(assays=list(counts=SRX_assay),
                         rowData=rowData(x),
                         colData=SRX_coldata,
                         metadata=metadata(x)))
}
