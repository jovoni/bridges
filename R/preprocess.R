
process_input_data = function(data,
                              chromosomes = c(1:22, "X", "Y"),
                              alleles = c("CN", "A", "B"),
                              k_jitter_fix = 2,
                              fillna = 0) {

  chromosomes = chromosomes[chromosomes %in% data$chr]
  alleles = alleles[alleles %in% colnames(data)]

  all_input_Xs = lapply(chromosomes, function(chromosome) {
    #print(chromosome)
    input_Xs = lapply(alleles, function(allele) {
      X = tibble_to_matrix(data, chromosome = chromosome, value_column = allele)
      X[is.na(X)] = fillna
      input_X = smooth_input_X(X, k_jitter_fix = k_jitter_fix)
      input_X
    })
    names(input_Xs) = alleles
    input_Xs
  })
  names(all_input_Xs) = chromosomes
  all_input_Xs
}

tibble_to_matrix <- function(data, chromosome, value_column = "CN") {
  # Filter data by chromosome
  data <- data %>%
    dplyr::filter(.data$chr == chromosome) %>%
    dplyr::select(.data$cell_id, .data$chr, .data$start, .data$end, value_column)

  # Check if value_column is valid
  if (!value_column %in% c("A", "B", "CN")) {
    stop("value_column must be one of 'A', 'B', or 'CN'")
  }

  # Get unique cell_ids for rows (maintain original order)
  unique_cells <- unique(data$cell_id)

  # Create a data frame with all unique (start, end) pairs for columns
  bins <- data[, c("start", "end")] %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$start, .data$end)

  # Create column names in the same format as the original function
  bin_names <- paste0(bins$start, "-", bins$end)

  # Initialize the matrix with NAs
  result_matrix <- matrix(
    NA,
    nrow = length(unique_cells),
    ncol = nrow(bins),
    dimnames = list(
      unique_cells,
      bin_names
    )
  )

  # Create indices for direct matrix assignment (much faster than looping)
  cell_indices <- match(data$cell_id, unique_cells)
  bin_indices <- match(paste0(data$start, "-", data$end), bin_names)

  # Create an index matrix where each row identifies a position in the result matrix
  idx <- cbind(cell_indices, bin_indices)

  # Fill the matrix with values in one operation
  result_matrix[idx] <- data[[value_column]]

  return(result_matrix)
}

smooth_input_X = function(X, k_jitter_fix) {
  if (k_jitter_fix != 0) {
    X = jitter_fix(X, k = k_jitter_fix)
  }
  input_X = reduce_X(X, find_changing_columns(X) - 1)
  input_X
}

jitter_fix = function(X, k=2) {
  y = apply(X, MARGIN = 1, function(row) {
    as.numeric(abs(diff(row)) != 0)
  }) %>% t()
  y_column_names = colnames(y) = colnames(X)[2:ncol(X)]

  column_queue = colMeans(y) %>% sort(decreasing = TRUE)
  column_visited = c()

  c = names(column_queue)[1]
  for (c in names(column_queue)) {
    c_idx = which(y_column_names == c)
    neighbours_indexes = (c_idx-k):(c_idx+k)
    neighbours_indexes = neighbours_indexes[neighbours_indexes > 0 & neighbours_indexes != c_idx & neighbours_indexes <= ncol(y)]
    for (n_idx in neighbours_indexes) {
      if (!n_idx %in% column_visited) {
        y[,c_idx] = y[,c_idx] | y[,n_idx]
        y[,n_idx] = 0
        column_visited = c(column_visited, n_idx)
      }
    }
  }

  # Start with a copy of X
  X_smoothed <- X

  # Process each row
  for (i in 1:nrow(X)) {
    # Get breakpoint positions for this row
    breakpoints <- which(y[i,] != 0)

    # Add start and end points to create complete segments
    all_points <- c(1, breakpoints, ncol(X) + 1)

    # Pre-allocate a vector for the smoothed values in this row
    smoothed_row <- X[i,]

    # Calculate segment averages in one step
    for (j in 1:(length(all_points) - 1)) {
      start_idx <- all_points[j]
      end_idx <- all_points[j+1] - 1

      # Calculate segment average and assign to the entire segment at once
      smoothed_row[start_idx:end_idx] <- mean(X[i, start_idx:end_idx])
    }

    # Assign the entire smoothed row at once
    X_smoothed[i,] <- smoothed_row
  }

  X_smoothed = round(X_smoothed)
  X_smoothed
}

find_changing_columns <- function(matrix_data) {
  # Convert to matrix if not already
  if (!is.matrix(matrix_data)) {
    matrix_data <- as.matrix(matrix_data)
  }

  # Check if matrix contains integers
  if (!all(matrix_data == floor(matrix_data))) {
    warning("Matrix contains non-integer values. Converting to integers.")
    matrix_data <- floor(matrix_data)
  }

  # Number of columns
  ncols <- ncol(matrix_data)

  # If there's only one column, no changes to detect
  if (ncols <= 1) {
    return(integer(0))
  }

  # Initialize result vector
  change_columns <- integer(0)

  # Check each pair of adjacent columns
  for (j in 2:ncols) {
    # Get the current and previous columns
    current_col <- matrix_data[, j]
    prev_col <- matrix_data[, j-1]

    # Check if any value changes between columns
    if (any(current_col != prev_col)) {
      # Store the column index where the change was detected
      change_columns <- c(change_columns, j)
    }
  }

  return(change_columns)
}

reduce_X = function(X, breakpoints) {
  sorted_bps = sort(breakpoints)
  sorted_bps = c(0, sorted_bps, ncol(X))

  original_colnames = colnames(X)

  reduce_row = function(row) {
    lapply(2:length(sorted_bps), function(i) {
      idxs = (sorted_bps[i-1]+1):(sorted_bps[i])
      round(mean(row[idxs]))
    }) %>% unlist()
  }

  new_X = lapply(1:nrow(X), function(j){
    reduce_row(X[j,])
  }) %>% do.call("rbind", .)
  rownames(new_X) = rownames(X)

  # Now assign new column names based on segment ranges
  new_colnames = sapply(2:length(sorted_bps), function(i) {
    start_idx = sorted_bps[i-1] + 1
    end_idx = sorted_bps[i]
    start_pos = sub("-.*", "", original_colnames[start_idx])
    end_pos = sub(".*-", "", original_colnames[end_idx])
    paste0(start_pos, ":", end_pos)
  })

  colnames(new_X) = new_colnames
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
