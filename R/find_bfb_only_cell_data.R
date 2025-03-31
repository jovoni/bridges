find_bfb_only_cell_data <- function(cell_data) {
  # Add root if not present
  if (!"root" %in% cell_data$cell_id) {
    root_row <- dplyr::tibble(
      cell_id = "root",
      parent_id = NA,
      bfb_event = FALSE,
      birth_time = 0,
      death_time = 0
    )
    cell_data <- dplyr::bind_rows(root_row, cell_data)
  }

  cell_data$parent_id[is.na(cell_data$parent_id) & cell_data$cell_id != "root"] <- "root"

  # Function to find the most recent BFB progenitor
  find_bfb_progenitor <- function(cell_id) {
    current <- cell_id
    while (!is.na(current)) {
      current_row = cell_data %>% dplyr::filter(cell_id == current)
      if (current_row$bfb_event) return(current_row$parent_id)
      parent_row = cell_data %>% dplyr::filter(cell_id == current_row$parent_id)
      if (nrow(parent_row) == 0 || is.na(parent_row$parent_id)) return(NA)
      if (parent_row$bfb_event) return(parent_row$cell_id)
      current <- current_row$parent_id
    }
    return("root")  # Default to root if no BFB ancestor found
  }

  # Update parent_id to reflect the most recent BFB progenitor
  cell_data = cell_data %>%
    dplyr::mutate(parent_id = sapply(cell_id, find_bfb_progenitor))

  cell_data %>%
    dplyr::filter(is_alive | cell_id %in% cell_data$parent_id)
}
