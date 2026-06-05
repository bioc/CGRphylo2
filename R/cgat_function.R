# Suppress NOTE about global variable used by cgrplot() and plot_cgr()
utils::globalVariables(c("fasta_filtered"))

# Internal version used by parallelCGR — takes sequences directly to avoid
# writing to the global environment.
cgat_local <- function(k_mer, seq_index, len_trim, sequences) {
  sequence <- sequences[[seq_index]]
  sequence <- substr(sequence, 1, len_trim)
  sequence <- toupper(sequence)
  sequence <- gsub("[^ACGT]", "", sequence)
  n_kmers <- 4^k_mer
  bases <- c("A", "C", "G", "T")
  all_kmers <- apply(expand.grid(rep(list(bases), k_mer)), 1, paste, collapse = "")
  freq_matrix <- matrix(0, nrow = n_kmers, ncol = 1)
  rownames(freq_matrix) <- all_kmers
  seq_length <- nchar(sequence)
  if (seq_length >= k_mer) {
    for (i in seq_len(seq_length - k_mer + 1)) {
      kmer <- substr(sequence, i, i + k_mer - 1)
      if (kmer %in% all_kmers) {
        freq_matrix[kmer, 1] <- freq_matrix[kmer, 1] + 1
      }
    }
    total_kmers <- sum(freq_matrix)
    if (total_kmers > 0) freq_matrix <- freq_matrix / total_kmers
  }
  return(freq_matrix)
}


#' Calculate distance between two CGR frequency matrices
#'
#' Computes the distance between two CGR frequency matrices using the specified
#' distance metric.
#'
#' @param matrix1 Numeric matrix. First CGR frequency matrix.
#' @param matrix2 Numeric matrix. Second CGR frequency matrix.
#' @param distance_type Character. Type of distance to calculate. Options are:
#'   \itemize{
#'     \item "Euclidean" (default): Standard Euclidean distance
#'     \item "S_Euclidean": Squared Euclidean distance
#'     \item "Manhattan": Manhattan (city block) distance
#'   }
#'
#' @return Numeric. The calculated distance between the two matrices.
#'
#' @details
#' This function calculates pairwise distances between CGR frequency matrices.
#' The Euclidean distance is most commonly used, but Manhattan and squared
#' Euclidean distances are also available for specific applications.
#'
#' @examples
#' # Create two simple frequency matrices
#' matrix1 <- matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 1)
#' matrix2 <- matrix(c(0.15, 0.25, 0.25, 0.35), ncol = 1)
#'
#' # Calculate Euclidean distance
#' dist_euclidean <- matrixDistance(matrix1, matrix2,
#'                                  distance_type = "Euclidean")
#' print(dist_euclidean)
#'
#' # Calculate Manhattan distance
#' dist_manhattan <- matrixDistance(matrix1, matrix2,
#'                                  distance_type = "Manhattan")
#' print(dist_manhattan)
#'
#' @export
matrixDistance <- function(matrix1, matrix2, distance_type = "Euclidean") {
  if (distance_type == "Euclidean") {
    dist <- sqrt(sum((matrix1 - matrix2)^2))
  } else if (distance_type == "S_Euclidean") {
    dist <- sum((matrix1 - matrix2)^2)
  } else if (distance_type == "Manhattan") {
    dist <- sum(abs(matrix1 - matrix2))
  } else {
    stop("Invalid distance_type. Choose 'Euclidean', 'S_Euclidean', or 'Manhattan'")
  }
  return(dist)
}


#' Filter sequences by N content
#'
#' Filters a list of DNA sequences by removing those with too many ambiguous
#' (N) bases.
#'
#' @param fastafile List. A named list of DNA sequences (e.g. from
#'   \code{seqinr::read.fasta()}).
#' @param N_filter Integer. Maximum number of N bases allowed in a sequence.
#'   Sequences with more N's than this threshold will be removed.
#'
#' @return List. Filtered list of sequences.
#'
#' @details
#' This function is useful for quality control before phylogenetic analysis.
#' Sequences with excessive ambiguous bases can affect the accuracy of distance
#' calculations and tree construction.
#'
#' @examples
#' test_seqs <- list(
#'     good_seq = "ATCGATCG",
#'     bad_seq  = "ATCGNNNNNATCG",
#'     okay_seq = "ATCGNNATCG"
#' )
#' filtered <- filter_N(test_seqs, N_filter = 3)
#' length(filtered)  # Should be 2
#'
#' @export
filter_N <- function(fastafile, N_filter) {
  n_counts <- vapply(fastafile, function(seq) {
    length(grep("N", strsplit(toupper(seq), "")[[1]]))
  }, FUN.VALUE = integer(1))
  filtered <- fastafile[n_counts <= N_filter]
  message(
    "Filtered ", length(fastafile) - length(filtered),
    " sequences with >", N_filter, " N bases"
  )
  return(filtered)
}

#' Create metadata table for sequences
#'
#' Extracts and compiles metadata including sequence length, GC content, and
#' N content from FASTA sequences.
#'
#' @param fastafile List. A list of DNA sequences read by seqinr::read.fasta()
#' @param N_filter Integer. N filter threshold (for reference in output)
#'
#' @return data.frame. A data frame with columns: name, length, GC_content,
#'   N_content
#'
#' @details
#' This function provides useful summary statistics for quality control and
#' understanding sequence characteristics before phylogenetic analysis.
#'
#' @examples
#' # Create test sequences
#' test_seqs <- list(
#'     seq1 = "ATCGATCGATCG",
#'     seq2 = "GCTAGCTAGCTA",
#'     seq3 = "AAAATTTTCCCCGGGG"
#' )
#'
#' # Get metadata
#' meta <- create_meta(test_seqs, N_filter = 50)
#' print(meta)
#'
#' @export
create_meta <- function(fastafile, N_filter) {
  meta_df <- data.frame(
    name = names(fastafile),
    length = vapply(fastafile, nchar, FUN.VALUE = integer(1)),
    GC_content = vapply(fastafile, function(seq) {
      bases <- strsplit(toupper(seq), "")[[1]]
      gc_count <- sum(bases %in% c("G", "C"))
      gc_count / length(bases) * 100
    }, FUN.VALUE = numeric(1)),
    N_content = vapply(fastafile, function(seq) {
      bases <- strsplit(toupper(seq), "")[[1]]
      as.integer(sum(bases == "N"))
    }, FUN.VALUE = integer(1)),
    stringsAsFactors = FALSE
  )
  return(meta_df)
}
