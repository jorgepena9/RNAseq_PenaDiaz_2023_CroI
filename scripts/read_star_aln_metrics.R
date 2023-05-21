###------------------------------------------------------------------------

#' Import STAR Log.final.out files for analyzing mapping rates.
#'
#' \code{read_star_aln_metrics} returns a dataframe from STAR Log.final.out files.
#'
#' This function will read in all the '*Log.final.out' files in the specified
#' directory and returns a dataframe with alignment metrics.
#'
#' @param dir A character scalar indicating the directory containing all the
#'     STAR Log.final.out files.
#' @param rgx A character scalar representing a regex used to parse out the
#'     sample name from the name of the Log.final.out files.
#' @return A dataframe.
#' @importFrom dplyr mutate across
#' @importFrom tibble column_to_rownames
#' @importFrom tidyselect ends_with
#' @export
read_star_aln_metrics <- function(dir, rgx){
  
  # get all the file names for the aln metrics files from STAR
  aln_metrics_files <- grep(".*Log.final.out(\\.[^\\.\\s]+)?$",
                            list.files(dir, full.names = TRUE),
                            value = TRUE,
                            perl = TRUE)
  
  stopifnot(length(aln_metrics_files) > 0) # check at least 1 file
  
  # read in files
  metrics_of_interest <- c(input_reads =
                             "Number of input reads",
                           uniq_aln_reads =
                             "Uniquely mapped reads number",
                           mult_aln_reads =
                             "Number of reads mapped to multiple loci")
  
  # parse out rows of interest
  aln_metrics_list <- lapply(aln_metrics_files,
                             function(x){
                               samp_name <- stringr::str_extract(basename(x),
                                                                 rgx)
                               
                               metrics <- readr::read_delim(x,
                                                            delim = "\\n",
                                                            col_names = FALSE,
                                                            col_types = "c")
                               
                               stopifnot(ncol(metrics) == 1)
                               
                               keep_rows_list <- lapply(metrics_of_interest,
                                                        grep,
                                                        x = metrics[[1]])
                               
                               keep_rows <- metrics[unlist(keep_rows_list), ]
                               
                               keep_metrics <- stringr::str_extract(
                                 keep_rows[[1]],
                                 "(?<=\\t)\\d+$")
                               
                               names(keep_metrics) <- names(metrics_of_interest)
                               
                               keep_metrics["sample"] <- samp_name
                               
                               keep_metrics_df <- data.frame(
                                 as.list(keep_metrics),
                                 stringsAsFactors = FALSE)
                               
                               return(keep_metrics_df)
                             })
  
  aln_metrics_mat <- do.call(rbind, aln_metrics_list)
  
  aln_metrics_mat <- aln_metrics_mat %>%
    dplyr::mutate(dplyr::across(tidyselect::ends_with("_reads"), as.numeric)) %>%
    tibble::column_to_rownames(var="sample")
  
  aln_metrics_mat <- data.matrix(aln_metrics_mat)
  
  return(aln_metrics_mat)
}