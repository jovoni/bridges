
process_input_data <- function(data,
                                    chromosomes = c(1:22, "X", "Y"),
                                    alleles = c("CN", "A", "B"),
                                    k_jitter_fix = 2,
                                    fillna = 0) {
  chromosomes <- intersect(chromosomes, unique(data$chr))
  alleles     <- intersect(alleles, colnames(data))

  # Split once by chromosome to cut repeated filtering
  split_by_chr <- split(data, data$chr)
  split_by_chr <- split_by_chr[intersect(names(split_by_chr), as.character(chromosomes))]

  worker <- function(df_chr) {
    out <- lapply(alleles, function(allele) {
      X <- tibble_to_matrix(df_chr, value_column = allele)
      X[is.na(X)] <- fillna
      smooth_input_X(X, k_jitter_fix = k_jitter_fix)
    })
    names(out) <- alleles
    out
  }

  all_input_Xs <- lapply(split_by_chr, worker)
  # Keep original order of chromosomes requested
  all_input_Xs <- all_input_Xs[as.character(chromosomes)]
  names(all_input_Xs) <- as.character(chromosomes)
  all_input_Xs
}


tibble_to_matrix <- function(data_chr, value_column = "CN") {
  if (!value_column %in% c("A", "B", "CN")) {
    stop("value_column must be one of 'A', 'B', or 'CN'")
  }

  # Use data.table for much faster operations
  if (!inherits(data_chr, "data.table")) {
    data_chr <- data.table::as.data.table(data_chr)
  }

  # Get unique values more efficiently
  unique_cells <- data_chr[, unique(cell_id)]
  bins_dt <- data_chr[, unique(.SD), .SDcols = c("start", "end")]
  data.table::setorder(bins_dt, start, end)

  bin_names <- paste0(bins_dt$start, "-", bins_dt$end)

  # Pre-allocate matrix
  res <- matrix(
    NA_real_,
    nrow = length(unique_cells),
    ncol = nrow(bins_dt),
    dimnames = list(unique_cells, bin_names)
  )

  # Use data.table for faster matching and assignment
  data_chr[, bin_name := paste0(start, "-", end)]
  data_chr[, row_idx := match(cell_id, unique_cells)]
  data_chr[, col_idx := match(bin_name, bin_names)]

  # Vectorized assignment - much faster than cbind indexing
  valid_rows <- !is.na(data_chr$row_idx) & !is.na(data_chr$col_idx)
  if (any(valid_rows)) {
    indices <- (data_chr$col_idx[valid_rows] - 1L) * length(unique_cells) +
      data_chr$row_idx[valid_rows]
    res[indices] <- data_chr[[value_column]][valid_rows]
  }

  res
}


smooth_input_X <- function(X, k_jitter_fix) {
  if (k_jitter_fix != 0) {
    X <- jitter_fix(X, k = k_jitter_fix)
  }
  reduce_X(X, breakpoints = find_changing_columns(X) - 1L)
}

jitter_fix <- function(X, k = 2) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n_rows <- nrow(X)
  n_cols <- ncol(X)

  if (n_cols <= 1) return(round(X))

  # Step 1: Vectorized difference calculation (faster than apply)
  y <- (abs(X[, -1, drop = FALSE] - X[, -n_cols, drop = FALSE]) != 0) * 1.0

  # Set column names to match original
  y_column_names <- colnames(X)[2:n_cols]
  colnames(y) <- y_column_names

  # Step 2: Create priority queue - same as original
  column_queue <- sort(colMeans(y), decreasing = TRUE)

  # Step 3: Optimized neighbor processing but same logic
  # Pre-allocate visited vector instead of growing it
  column_visited <- logical(ncol(y))
  names(column_visited) <- 1:ncol(y)

  for (c in names(column_queue)) {
    c_idx <- which(y_column_names == c)

    # Same neighbor calculation as original
    neighbours_indexes <- (c_idx - k):(c_idx + k)
    neighbours_indexes <- neighbours_indexes[neighbours_indexes > 0 &
                                               neighbours_indexes != c_idx &
                                               neighbours_indexes <= ncol(y)]

    # Process neighbors with same logic but faster indexing
    for (n_idx in neighbours_indexes) {
      if (!column_visited[n_idx]) {
        # Same OR operation and zeroing as original
        y[, c_idx] <- y[, c_idx] | y[, n_idx]
        y[, n_idx] <- 0
        column_visited[n_idx] <- TRUE
      }
    }
  }

  # Step 4: Optimized smoothing but same algorithm
  X_smoothed <- X

  # Pre-compute cumulative sums for faster segment means
  for (i in 1:nrow(X)) {
    breakpoints <- which(y[i, ] != 0)

    if (length(breakpoints) == 0) {
      # No breakpoints - smooth entire row
      X_smoothed[i, ] <- mean(X[i, ])
    } else {
      # Same segment logic as original
      all_points <- c(1, breakpoints, ncol(X) + 1)

      # Use cumulative sum for faster mean calculation
      row_cumsum <- c(0, cumsum(X[i, ]))

      for (j in 1:(length(all_points) - 1)) {
        start_idx <- all_points[j]
        end_idx <- all_points[j + 1] - 1

        # Faster mean calculation using cumsum
        segment_mean <- (row_cumsum[end_idx + 1] - row_cumsum[start_idx]) /
          (end_idx - start_idx + 1)
        X_smoothed[i, start_idx:end_idx] <- segment_mean
      }
    }
  }

  round(X_smoothed)
}


find_changing_columns <- function(X) {
  if (!is.matrix(X)) X <- as.matrix(X)
  # If integer-like not guaranteed, skip the expensive check and just use direct comparison
  if (ncol(X) <= 1L) return(integer(0))
  # Compare adjacent columns in one shot
  changes <- X[, -1, drop = FALSE] != X[, -ncol(X), drop = FALSE]
  which(colSums(changes) > 0) + 1L
}

reduce_X = function(X, breakpoints) {
  sorted_bps = sort(breakpoints)
  sorted_bps = c(0, sorted_bps, ncol(X))
  original_colnames = colnames(X)

  # Create segment matrix more efficiently
  n_segments = length(sorted_bps) - 1

  # Use sapply for vectorized segment processing
  new_X = sapply(1:n_segments, function(i) {
    cols = (sorted_bps[i] + 1):sorted_bps[i + 1]
    if (length(cols) == 1) {
      round(X[, cols])
    } else {
      round(rowMeans(X[, cols, drop = FALSE]))
    }
  })

  # Ensure it's a matrix
  if (!is.matrix(new_X)) {
    new_X = as.matrix(new_X)
  }

  rownames(new_X) = rownames(X)

  # Vectorized column names
  colnames(new_X) = sapply(1:n_segments, function(i) {
    start_idx = sorted_bps[i] + 1
    end_idx = sorted_bps[i + 1]
    start_pos = sub("-.*", "", original_colnames[start_idx])
    end_pos = sub(".*-", "", original_colnames[end_idx])
    paste0(start_pos, ":", end_pos)
  })

  new_X
}


create_diploid_data = function(data) {
  cid = data$cell_id[1] # sample any cell_id
  data %>%
    dplyr::filter(cell_id == cid) %>%
    dplyr::mutate(
      cell_id = "diploid",
      A = ifelse(chr %in% c("X", "Y"), 1, 1),
      B = ifelse(chr %in% c("X", "Y"), 0, 1),
      CN = A + B
    )
}
