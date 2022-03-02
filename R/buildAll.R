
# hsapiens column data we are targetting; notably the accession are in cols$SRR_accession
cols <- DataFrame(read.csv(system.file("inst", "hsapiens_colData_transitions_v3.5.csv", package="homosapienDEE2CellScore")))

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

buildRaw <- function(species="hsapiens", accessions=as.list(unique(cols$SRR_accession)), quiet=TRUE, metadata=getDEE2Metadata(species, quiet=quiet)) {
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

buildData <- function(species="hsapiens", name_prefix="homosapienDEE2Data", name_suffix=".csv", build_raw=FALSE, build_srx_agg=FALSE, build_deseq2=TRUE, build_tsne=TRUE, build_rank=TRUE, generate_qc_pass = TRUE, generate_qc_warn = TRUE, base=getwd(), quiet=TRUE, metadata=if(!(build_raw || build_srx_agg || build_deseq2 || build_tsne || build_rank) || !(generate_qc_pass || generate_qc_warn)) { return(list()); } else { getDEE2Metadata(species, quiet=quiet) }, counts.cutoff = 10, accessions=as.list(unique(cols$SRR_accession)), in_data = if(!(build_raw || build_srx_agg || build_deseq2 || build_tsne || build_rank) || !(generate_qc_pass || generate_qc_warn)) { return(list()); } else { buildRaw(species=species, accessions=accessions, quiet=quiet, metadata=metadata) }, dds_design = ~ 1, write_files = TRUE) {

  out <- list()
  outputs <- list()
  # Check whether we are not going to do something
  if(!(build_raw || build_srx_agg || build_deseq2 || build_tsne) || !(generate_qc_pass || generate_qc_warn)) {
    return(out);
  }
  # All of the 'optionally overriden' data possible is calculated in the function's arguments.

  # First we add on the colData from our pre-curated package colData
  with_colData <- in_data
  print(paste("with_colData nrows:", nrow(with_colData), "ncols:", ncol(with_colData)))
  colData_cols <- cols
  rownames(colData_cols) <- colData_cols$SRR_accession
  print(paste("colData_cols nrows:", nrow(colData_cols), "ncols:", ncol(colData_cols)))
  colData_cols <- colData_cols[colData_cols$SRR_accession %in% rownames(colData(with_colData)),]
  colData_cols <- colData_cols[!duplicated(colData_cols$SRR_accession),]
  print(paste("deduplicated colData_cols nrows:", nrow(colData_cols), "ncols:", ncol(colData_cols)))
  colData(with_colData)<-colData_cols  #cols[!duplicated(cols$SRR_accession),]

  # Take either the clean pass of qc data, or the data which passes and the data which has warnings, but not failures
  qc_pass <- with_colData[, startsWith(with_colData$QC_summary, "PASS")]
  qc_warn <- with_colData[, startsWith(with_colData$QC_summary, "PASS") | startsWith(with_colData$QC_summary, "WARN")]

  if(build_raw) {
    if(generate_qc_pass) {
      out <- c(out, list(qc_pass_raw=.finishUp(qc_pass)))
      outputs <- c(outputs, list(qc_pass_raw=paste(name_prefix, "_PASS_raw", sep="")))
    }
    if(generate_qc_warn) {
      out <- c(out, list(qc_warn_raw=.finishUp(qc_warn)))
      outputs <- c(outputs, list(qc_warn_raw=paste(name_prefix, "_WARN_raw", sep="")))
    }
  }

  # Now we aggregate into srx, instead of by srr
  qc_pass_agg <- srx_agg_se(qc_pass)
  qc_warn_agg <- srx_agg_se(qc_warn)

  # Now we filter based on gene activity
  qc_pass_filtered <- qc_pass_agg[rowSums(assay(qc_pass_agg)) > counts.cutoff,]
  qc_warn_filtered <- qc_warn_agg[rowSums(assay(qc_warn_agg)) > counts.cutoff,]
  calls_pass_template <- assay(addCallData(qc_pass_filtered), "calls")
  calls_warn_template <- assay(addCallData(qc_warn_filtered), "calls")

  if(build_srx_agg) {
    if(generate_qc_pass) {
      out <- c(out, list(qc_pass_agg=.finishUp(qc_pass_filtered)))
      outputs <- c(outputs, list(qc_pass_agg=paste(name_prefix, "_PASS_agg", sep="")))
    }
    if(generate_qc_warn) {
      out <- c(out, list(qc_warn_agg=.finishUp(qc_warn_filtered)))
      outputs <- c(outputs, list(qc_warn_agg=paste(name_prefix, "_WARN_agg", sep="")))
    }
  }

  if(build_deseq2) {
    # Now normalisation
    if(generate_qc_pass) {
      dds_qc_pass_filtered <- BiocGenerics::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(qc_pass_filtered, "counts"),
        colData = SummarizedExperiment::colData(qc_pass_filtered),
        rowData = SummarizedExperiment::rowData(qc_pass_filtered),
        design = dds_design))
      logcounts_qc_pass_filtered <- log2(counts(dds_qc_pass_filtered, normalize=TRUE) + 1)
      deseq2_final <- SummarizedExperiment(assays=list(counts=logcounts_qc_pass_filtered),
                                          rowData=SummarizedExperiment::rowData(qc_pass_filtered),
                                          colData=SummarizedExperiment::colData(qc_pass_filtered),
                                          metadata=metadata(qc_pass_filtered))
      out <- c(out, list(qc_pass_deseq2=.finishUp(deseq2_final)))
      outputs <- c(outputs, list(qc_pass_deseq2=paste(name_prefix, "_PASS_deseq2", sep="")))
    }
    if(generate_qc_warn) {
      dds_qc_warn_filtered <- BiocGenerics::estimateSizeFactors(DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(qc_warn_filtered, "counts"),
        colData = SummarizedExperiment::colData(qc_warn_filtered),
        rowData = SummarizedExperiment::rowData(qc_warn_filtered),
        design = dds_design))
      logcounts_qc_warn_filtered <- log2(counts(dds_qc_warn_filtered, normalize=TRUE) + 1)
      deseq2_final <- SummarizedExperiment(assays=list(counts=logcounts_qc_warn_filtered),
                                          rowData=SummarizedExperiment::rowData(qc_warn_filtered),
                                          colData=SummarizedExperiment::colData(qc_warn_filtered),
                                          metadata=metadata(qc_warn_filtered))
      out <- c(out, list(qc_warn_deseq2=.finishUp(deseq2_final)))
      outputs <- c(outputs, list(qc_warn_deseq2=paste(name_prefix, "_WARN_deseq2", sep="")))
    }
  }
  # We don't add calls data to the tsne form yet; I need to check whether that will work properly
  if(build_tsne) {
    if(generate_qc_pass) {
      tsne_qc_pass_filtered <- Rtsne(unique(t(as.matrix(SummarizedExperiment::as.data.frame(SummarizedExperiment::assay(qc_pass_filtered, "counts"))))), check_duplicates=FALSE)$Y
      out <- c(out, list(qc_pass_tsne=tsne_qc_pass_filtered))
      outputs <- c(outputs, list(qc_pass_tsne=paste(name_prefix, "_PASS_tsne", sep="")))
    }
    if(generate_qc_warn) {
      tsne_qc_warn_filtered <- Rtsne(unique(t(as.matrix(SummarizedExperiment::as.data.frame(SummarizedExperiment::assay(qc_warn_filtered, "counts"))))), check_duplicates=FALSE)$Y
      out <- c(out, list(qc_warn_tsne=tsne_qc_warn_filtered))
      outputs <- c(outputs, list(qc_warn_tsne=paste(name_prefix, "_WARN_tsne", sep="")))
    }
  }
  if (build_rank) {
    if(generate_qc_pass) {
      len <- length(assay(qc_pass_filtered, "counts"))
      #ranks <- lapply(DataFrame(t(assay(qc_pass_filtered, "counts"))), function (x) { (order(x, decreasing=TRUE)/len)})
      #ranks <- t(as.matrix((DataFrame(lapply(DataFrame(t(assay(qc_pass_filtered, "counts")[,])), function (x) { (order(x, decreasing=TRUE)/length(x))}), row.names=colnames(assay(qc_pass_filtered, "counts")[,])))))
      nrows <- length(row.names(assay(qc_pass_filtered, "counts")))
      #ranks <- (nrows - rowRanks(assay(qc_pass_filtered, "counts"), ties.method="average") + 1) / nrows
      ranks <- colRanks(assay(qc_pass_filtered, "counts"), ties.method="average", preserveShape=TRUE) / nrows
      rank <- SummarizedExperiment(assays=list(counts=ranks, calls=calls_pass_template),
                                          rowData=SummarizedExperiment::rowData(qc_pass_filtered),
                                          colData=SummarizedExperiment::colData(qc_pass_filtered),
                                          metadata=metadata(qc_pass_filtered))
      out <- c(out, list(qc_pass_rank=addProbeId(rank)))
      outputs <- c(outputs, list(qc_pass_rank=paste(name_prefix, "_PASS_rank", sep="")))
    }
    if(generate_qc_warn) {
      len <- length(assay(qc_warn_filtered, "counts"))
      #ranks <- lapply(DataFrame(t(assay(qc_warn_filtered, "counts"))), function (x) { (order(x, decreasing=TRUE)/len)})
      #ranks <- t(as.matrix((DataFrame(lapply(DataFrame(t(assay(qc_warn_filtered, "counts")[,])), function (x) { (order(x, decreasing=TRUE)/length(x))}), row.names=colnames(assay(qc_warn_filtered, "counts")[,])))))
      nrows <- length(row.names(assay(qc_warn_filtered, "counts")))
      #ranks <- (nrows - rowRanks(assay(qc_warn_filtered, "counts"), ties.method="average") + 1) / nrows
      ranks <- colRanks(assay(qc_warn_filtered, "counts"), ties.method="average", preserveShape=TRUE) / nrows
      rank <- SummarizedExperiment(assays=list(counts=ranks, calls=calls_warn_template),
                                          rowData=SummarizedExperiment::rowData(qc_warn_filtered),
                                          colData=SummarizedExperiment::colData(qc_warn_filtered),
                                          metadata=metadata(qc_warn_filtered))
      out <- c(out, list(qc_warn_rank=addProbeId(rank)))
      outputs <- c(outputs, list(qc_warn_rank=paste(name_prefix, "_WARN_rank", sep="")))
    }
  }

  if (write_files) {
    lapply(names(out), function(n) {
      writeOutSE(the_summarized_experiment=out[[n]], filename_base=outputs[[n]], filename_ext=name_suffix)
    })
  }
  out
}

# Takes a named-list of processed data plus a bunch of filenames for each named thing
writeOutput <- function(the_data, outputs=list(dds_qc_pass_filtered="dds_qc_pass_filtered.csv", dds_qc_warn_filter="dds_qc_warn_filter.csv")) {
  lapply(names(outputs), function(n) {
    write.csv(the_data[[n]], file=outputs[[n]], row.names=TRUE)
  })
}

writeOutSE <- function(
                       the_summarized_experiment, filename_base="SE_out", filename_ext=".csv", filenames=list(
                       metadata=paste(filename_base, "_metadata", filename_ext,sep=""),
                       assay_counts=paste(filename_base, "_assay_counts", filename_ext,sep=""),
                       assay_calls=paste(filename_base, "_assay_calls", filename_ext,sep=""),
                       colData=paste(filename_base, "_colData", filename_ext,sep=""),
                       rowData=paste(filename_base, "_rowData", filename_ext,sep=""))
                      ) {
  the_metadata<-metadata(the_summarized_experiment)
  the_assay_counts<-assay(the_summarized_experiment, "counts")
  the_assay_calls<-assay(the_summarized_experiment, "calls")
  the_colData<-colData(the_summarized_experiment)
  the_rowData<-rowData(the_summarized_experiment)
  writeOutput(list(metadata=the_metadata, assay_counts=the_assay_counts, assay_calls=the_assay_calls, colData=the_colData, rowData=the_rowData), outputs=filenames)
}

readInSE <- function(metadata_file="SE_out_metadata.csv", assay_counts_file="SE_out_assay_counts.csv", assay_calls_file="SE_out_assay_calls.csv", colData_file="SE_out_colData.csv", rowData_file="SE_out_rowData.csv") {
  metadata_in <- list() #read.csv(metadata_file, row.names=TRUE)
  assay_counts_in <- read.csv(assay_counts_file, row.names=1)
  assay_calls_in <- read.csv(assay_calls_file, row.names=1)
  colData_in <- read.csv(colData_file)
  rowData_in <- read.csv(rowData_file, header=FALSE, row.names=1, fill=TRUE, skipNul=TRUE, blank.lines.skip=TRUE)
  #Hack to 'fix' rowData - this doesn't actually fix things, it just breaks them in a way I don't currently care about
  if ((length(colnames(rowData_in)) == 1) && (colnames(rowData_in)[[1]]=="V2")) {
    rowData_in<-rowData_in[NULL]
  }
  return(SummarizedExperiment(assays=list(counts=assay_counts_in, calls=assay_calls_in), rowData=rowData_in, colData=colData_in, metadata=metadata_in))
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
    SRX_cols1 <- sapply(X=srx_es, function(srx) {
        srrs<-rownames(mds)[which(mds[["SRX_accession"]] %in% srx)]
        #list(SummarizedExperiment::assay(x, "counts")[,srrs])
        paste(srrs, collapse=" ")
    })
    SRX_col1real <- DataFrame(matrix(SRX_cols1, dimnames=list(srx_es, list("Aggregated_From"))))
    SRX_cols2 <- cols[!duplicated(cols$SRX_accession),]
    SRX_cols2 <- SRX_cols2[SRX_cols2$SRX_accession %in% srx_es,]
    rownames(SRX_cols2) <- SRX_cols2$SRX_accession
    print(paste("c1 nrow:", nrow(SRX_col1real), "c2 nrow:", nrow(SRX_cols2)))
    SRX_coldata <- cbind(SRX_col1real, SRX_cols2)
    return(SummarizedExperiment(assays=list(counts=SRX_assay),
                         rowData=rowData(x),
                         colData=SRX_coldata,
                         metadata=metadata(x)))
}

# For RNASeq data we assume any non-zero value is a call
addCallData <- function(summarized_experiment, threshold=0) {
  calls <- ifelse(assay(summarized_experiment, "counts") > threshold, 1, 0)
  assay(summarized_experiment, "calls") <- calls
  return(summarized_experiment)
}
addProbeId <- function(summarized_experiment) {
  rd <- rowData(summarized_experiment)
  probe_id <- rownames(rd)
  feature_id <- probe_id
  # We add both probe_id and feature_id to deal with a case where CellScore switches to using the more general feature_id.
  # For our use though, they are identical.
  rowData(summarized_experiment) <- cbind(rd, probe_id, feature_id)
  return(summarized_experiment)
}
.finishUp <- function(summarized_experiment) {
  return(addCallData(addProbeId(summarized_experiment)))
}
