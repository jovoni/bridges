
keep_only_bfb_branches <- function(cell_history) {
  # Add root if not present
  if (!"root" %in% cell_history$cell_id) {
    root_row <- dplyr::tibble(
      cell_id = "root",
      parent_id = NA,
      bfb_event = FALSE,
      birth_time = 0,
      death_time = 0
    )
    cell_history <- dplyr::bind_rows(root_row, cell_history)
  }

  cell_history$parent_id[is.na(cell_history$parent_id) & cell_history$cell_id != "root"] <- "root"

  # Function to find the most recent BFB progenitor
  find_bfb_progenitor <- function(cell_id) {
    current <- cell_id
    while (!is.na(current)) {
      current_row = cell_history %>% dplyr::filter(.data$cell_id == current)
      if (current_row$bfb_event) return(current_row$parent_id)
      parent_row = cell_history %>% dplyr::filter(.data$cell_id == current_row$parent_id)
      if (nrow(parent_row) == 0 || is.na(parent_row$parent_id)) return(NA)
      if (parent_row$bfb_event) return(parent_row$cell_id)
      current <- current_row$parent_id
    }
    return("root")  # Default to root if no BFB ancestor found
  }

  # Update parent_id to reflect the most recent BFB progenitor
  cell_history = cell_history %>%
    dplyr::mutate(parent_id = sapply(.data$cell_id, find_bfb_progenitor))

  cell_history %>%
    dplyr::filter(.data$is_alive | .data$cell_id %in% cell_history$parent_id)
}
