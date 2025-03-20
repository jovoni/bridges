#' Process multiple sequences (cells) to get copy number runs
#' @param sequences A list of sequences, each represented as a list of segments
#' @param L Maximum position to consider
#' @return A tibble with copy number runs for all sequences
multi_seqs2CNruns <- function(sequences, L = NULL) {
  # Process each sequence and combine results
  result <- lapply(seq_along(sequences), function(seq_idx) {
    seq_runs <- seq2CNruns(sequences[[seq_idx]], L = L)
    seq_runs$cell_idx <- seq_idx
    seq_runs
  })

  # Combine all results into a single tibble
  dplyr::bind_rows(result)
}
