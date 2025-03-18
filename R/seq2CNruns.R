#' Create a tibble of copy number runs from a sequence represented as seq
#' @param seq A list of seq, each with start, end, and direction
#' @param L Maximum value to consider (defaults to max value in seq)
#' @return A tibble with columns segment_start, segment_end, and copy_number
seq2CNruns <- function(seq, L = NULL) {

  # Early exit for empty seq
  if (length(seq) == 0) {
    return(dplyr::tibble(segment_start = integer(0),
                  segment_end = integer(0),
                  copy_number = integer(0)))
  }

  # Find the maximum value if not provided
  if (is.null(L)) {
    L <- 0
    for (seg in seq) {
      L <- max(L, seg$start, seg$end)
    }
  }

  # Initialize copy number array (1 to L)
  cn_counts <- integer(L)

  # Count occurrences of each value in the sequence
  for (seg in seq) {
    start_val <- seg$start
    end_val <- seg$end
    direction <- seg$direction

    if (direction == 1) {
      # Increasing sequence
      for (val in start_val:end_val) {
        if (val >= 1 && val <= L) {
          cn_counts[val] <- cn_counts[val] + 1  # 1-based indexing
        }
      }
    } else if (direction == -1) {
      # Decreasing sequence
      for (val in start_val:end_val) {
        if (val >= 1 && val <= L) {
          cn_counts[val] <- cn_counts[val] + 1  # 1-based indexing
        }
      }
    } else {
      # Single value or constant sequence
      if (start_val >= 1 && start_val <= L) {
        cn_counts[start_val] <- cn_counts[start_val] + 1  # 1-based indexing
      }
    }
  }

  # Find runs of identical copy numbers
  result <- dplyr::tibble()

  if (length(cn_counts) > 0) {
    run_start <- 1  # Start from position 1 (1-based indexing)
    current_cn <- cn_counts[1]

    for (i in 2:length(cn_counts)) {
      if (cn_counts[i] != current_cn) {
        # End of a run - add to result
        result <- dplyr::bind_rows(result,
                            dplyr::tibble(segment_start = run_start,
                                   segment_end = i - 1,  # 1-based indexing
                                   copy_number = current_cn))
        run_start <- i  # 1-based indexing
        current_cn <- cn_counts[i]
      }
    }

    # Add the final run

    result <- dplyr::bind_rows(result,
                        dplyr::tibble(segment_start = run_start,
                               segment_end = L,  # Last position is L
                               copy_number = current_cn))
  }

  return(result)
}
