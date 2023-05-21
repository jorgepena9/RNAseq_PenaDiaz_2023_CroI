###### Star to mat #####

###------------------------------------------------------------------------

#' Convert STAR read counts to count matrix.
#'
#' \code{star_to_mat} returns a count matrix from STAR ReadsPerGene.out.tab
#' files.
#'
#' This function will read in all the '*ReadsPerGene.out.tab' files in the
#' specified directory and converts it to a counts matrix.
#'
#' @param dir A character scalar indicating the directory containing all the
#'     STAR ReadsPerGene.out.tab files.
#' @param rgx A character scalar representing a regex used to parse out the
#'     sample name from the name of the ReadsPerGene.out.tab file.
#' @param column An integer indicating the column to extract counts from.
#'     1 = unstranded, 2 = 1st read strand aligned with RNA, 3 = 2nd read
#'     strand aligned with RNA  with RNA
#' @param rm_ens_vers Logical indicating whether version number should be
#'     removed from gene ID (as indicated by trailing period and integers).
#' @return A matrix where the column names are the sample names and the row
#'     names are the gene names
#' @importFrom utils capture.output
#' @export

star_to_mat <- function(dir, rgx, column, rm_ens_vers = TRUE){
  
  # get all the file names for the count files from STAR
  count_file_rgx <- ".*ReadsPerGene.out.tab(\\.[^\\.\\s]+)?$"
  count_files <- grep(count_file_rgx,
                      list.files(dir, full.names = TRUE),
                      value = TRUE, perl = TRUE)
  
  # check at least 1 count file
  if (length(count_files) < 1){
    stop(paste0("No count files of regex ", count_file_rgx, "detected\n"))
  }
  
  # check that column requested is valid
  if (!column %in% c(1, 2, 3)){
    stop(paste0("Column must be between 1 and 3. Input was: ", column, "\n"))
  }
  
  # read in files
  counts_list <- lapply(count_files,
                        function(x){
                          column_names <- c("gene_id",
                                            "unstr",
                                            "str_r1",
                                            "str_r2")
                          
                          # read
                          df <- readr::read_tsv(file = x,
                                                col_names = column_names,
                                                col_types = "ciii")
                          
                          # add 'sample' column based on filename
                          samp_name <- stringr::str_extract(basename(x), rgx)
                          
                          df$sample = samp_name
                          
                          return(df)
                        })
  
  # check same columns in every file
  colnames_list <- lapply(counts_list, colnames)
  stopifnot(all(sapply(colnames_list[-1], identical, colnames_list[[1]])))
  
  # get the counts column indicated by function call
  counts_col_name <- colnames_list[[1]][1 + column]
  
  # process files
  counts_list <- lapply(counts_list, function(x){
    # guess correct strand selection
    nofeat_cts <- x[x$gene_id=="N_noFeature", , drop=FALSE]
    
    # heuristic method to guess strand
    if(nofeat_cts$str_r1 < nofeat_cts$str_r2){
      strand_guess <- 2
    } else{
      strand_guess <- 3
    }
    if(nofeat_cts$unstr < (0.5*nofeat_cts[, 1+strand_guess])){
      strand_guess <- 1
    }
    if(strand_guess != column){
      warning(paste0("Guessed strandedness was ", colnames(x)[1+strand_guess],
                     " (column ", strand_guess, ").\n"),
              paste0(utils::capture.output(print(as.data.frame(x[1:6, ]))),
                     collapse = "\n"))
    }
    
    
    df <- x %>%
      dplyr::filter(!.data$gene_id %in% c("N_unmapped",
                                          "N_multimapping",
                                          "N_noFeature",
                                          "N_ambiguous")) %>%
      dplyr::select(.data$gene_id,
                    dplyr::matches(counts_col_name),
                    .data$sample)
    
    # check for missing data
    stopifnot(all(complete.cases(df)))
    
    return(df)
    
  })
  
  # check each file has same # of genes
  samp_lens <- unlist(lapply(counts_list, function(x){
    num_row <- nrow(x)
    num_uniq_genes <- length(unique(x$gene_id))
    stopifnot(identical(num_row, num_uniq_genes)) # check gene names are unique
    
    return(num_row)
  }))
  samp_lens_uniq <- unique(unname(samp_lens))
  stopifnot(length(samp_lens_uniq) == 1 & samp_lens_uniq > 1)
  
  # rbind
  counts_rbind <- do.call(rbind, counts_list)
  
  # check each file has same genes
  ## split by gene
  split_by_genes <- split(counts_rbind, counts_rbind$gene_id)
  
  ## samples per gene
  split_by_genes_nrows <- unlist(lapply(split_by_genes, nrow))
  
  ## unique # of samples per gene
  split_by_genes_nrows_uniq <- unique(unname(split_by_genes_nrows))
  
  ## check that all genes have same number of samples and that # is equal
  ## to the number of input samples.
  stopifnot(length(split_by_genes_nrows_uniq) == 1 &
              identical(split_by_genes_nrows_uniq, length(count_files)))
  
  # spread to form count matrix
  count_mat <- tidyr::spread(data = counts_rbind,
                             key = .data$sample,
                             value = !!counts_col_name) %>%
    tibble::column_to_rownames(var = "gene_id")
  
  # check # of genes unchanged
  stopifnot(identical(nrow(count_mat), samp_lens_uniq))
  
  # remove version number from gene ids if rm_ens_vers
  if (rm_ens_vers) {
    new_gene_ids <- stringr::str_remove(rownames(count_mat), "\\.\\d+$")
    
    # check that removing the version # does not affect uniqueness.
    stopifnot(identical(length(unique(new_gene_ids)), length(rownames(count_mat))))
    
    # set new rownames
    rownames(count_mat) <- new_gene_ids
  }
  
  count_mat <- as.matrix(count_mat)
  
  return(count_mat)
}

