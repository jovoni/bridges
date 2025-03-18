#' Process multiple sequences (cells) to get copy number runs
#' @param sequences A list of sequences, each represented as a list of segments
#' @param L Maximum position to consider
#' @return A tibble with copy number runs for all sequences
multi_seqs2CNruns <- function(sequences, L = NULL) {

  # Process each sequence
  result <- lapply(1:length(sequences), function(seq_idx) {
    seq2CNruns(sequences[[seq_idx]], L = L) %>%
      dplyr::mutate(cell_idx = seq_idx)
  })

  do.call(dplyr::bind_rows, result)
}
