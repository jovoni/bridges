#' Reconstruct Phylogenetic Tree from Simulation Data
#'
#' This function reconstructs a phylogenetic tree from a simulation dataset by iteratively clustering cells
#' based on their genetic similarity. The process continues until a root cell containing all original cells is identified.
#'
#' @param sim A simulation object containing final cells and cell lifetimes.
#' @param L Lenght of "normal" sequence
#' @param alpha A numeric value controlling the weight of certain genetic features in the distance calculation (default: 0.1).
#' @param cn_weight A numeric weight applied to copy number variations in distance calculations (default: 1e-3).
#' @param eps A threshold for determining twin cells (default: 0.01).
#' @param remove_cells Logical, whether to remove redundant/twin cells (default: TRUE).
#' @param verbose Logical, whether to print iteration details (default: FALSE).
#'
#' @return A list containing:
#'   - `cell_df`: A dataframe with inferred cell lineages and events.
#'   - `pd_list`: A list mapping parents to their daughters in the phylogenetic tree.
#'
#' @export
reconstruct_tree = function(sim, L, alpha = .1, cn_weight = 1e-3, eps = .01, remove_cells = TRUE, verbose = FALSE) {
  # Extract final cells and cell names from simulation data
  cells = sim$final_cells
  cells_names = sim$cell_lifetimes %>% dplyr::filter(is_alive) %>% dplyr::pull(cell_id)
  original_cells = cells_names = lapply(cells_names, function(c) {paste0("|",c,"|")}) %>% unlist()
  removed_cells = TRUE
  N0 = length(cells)

  # Initialize cell matrix C
  C = cells2countvectors(cells, L)
  rownames(C) = cells_names

  # Initialize distance matrix D
  D = compute_distance_matrix(C, alpha = alpha, cn_weight = cn_weight)
  colnames(D) = rownames(D) = rownames(C)

  # Initialize mask matrix M (tracks which distances are valid)
  M = matrix(data = FALSE, nrow = N0, ncol = N0)
  diag(M) = TRUE
  colnames(M) = rownames(M) = rownames(C)

  # Initialize parent-daughters list
  pd_list = list()
  iter = 0
  found_root = FALSE

  while (!found_root) {
    if (verbose) print(iter)

    # Mask distances already used
    D[M] = Inf

    # Find smallest distance and group cells
    new_cells_and_daughters = NULL
    minimal_distances = unique(sort(D))
    min_dist_idx = 0
    while (is.null(new_cells_and_daughters)) {
      min_dist_idx = min_dist_idx + 1
      min_dist = minimal_distances[min_dist_idx]
      otus = propose_groups_of_cells(D, min_dist)
      new_cells_and_daughters = clean_group_of_cells(otus, C)
    }

    new_cells = new_cells_and_daughters$parents
    daughters = new_cells_and_daughters$daughters

    # Identify and remove twin cells
    if (remove_cells) {
      removable_cells_idx = lapply(1:nrow(new_cells), function(i) {
        if (all(compute_distance_matrix(C[daughters[[i]],], alpha = 1, cn_weight = 0) < eps)) {
          i
        }
      }) %>% unlist()
      removable_cells = lapply(removable_cells_idx, function(i) {daughters[[i]]}) %>% unlist() %>% unique()
    } else {
      removable_cells = NULL
    }

    # Update parent-daughter relationships
    names(daughters) = rownames(new_cells)
    pd_list = c(pd_list, daughters)

    # Update distance matrix with new cells
    D_new <- compute_distance_matrix_enhanced(new_cells, alpha = alpha, cn_weight = cn_weight)
    D_old_new <- compute_distance_matrix_enhanced(C, new_cells, alpha = alpha, cn_weight = cn_weight)
    D <- rbind(cbind(D, D_old_new), cbind(t(D_old_new), D_new))

    # Update cell matrix C
    C = rbind(C, new_cells)

    # Update mask matrix M
    M = extend_M(M = M, C = C)
    rownames(M) = colnames(M) = rownames(D) = colnames(D) = rownames(C)
    M = update_M(M, daughters, pd_list, all = FALSE)

    # Check if we found the root cell
    is_root = function(cell_name, N0) {
      cells_used = unlist(strsplit(gsub("\\|", "", cell_name), "-"))
      length(cells_used) == N0
    }

    removable_idx = rownames(C) %in% removable_cells
    removed_cells = rbind(removed_cells, C[removable_idx,])
    C = C[!removable_idx, ]
    M = M[!removable_idx, !removable_idx]
    D = D[!removable_idx, !removable_idx]

    found_root = any(unlist(lapply(rownames(new_cells), function(x) {is_root(x, N0)})))
    if (verbose) print(max(unlist(lapply(rownames(C), function(x) {length(unlist(strsplit(gsub("\\|", "", x), "-")))}))) / N0)

    iter = iter + 1
  }

  # Extract root node
  root = rownames(C)[unlist(lapply(rownames(C), function(x) {is_root(x, N0)}))]
  C = rbind(C, removed_cells)

  # Build phylogenetic tree
  phylogeny = build_phylogeny(pd_list, parent = root, parent_id = NA)

  # Add metadata to phylogeny
  phylogeny$type = ifelse(phylogeny$cell_id %in% original_cells, "observed", "pseudo")
  phylogeny$distance = lapply(1:nrow(phylogeny), function(i) {
    row = phylogeny[i,]
    if (is.na(row$parent_id)) return(NA)
    euclidean_distance(C[row$cell_id,], C[row$parent_id,])
  })
  phylogeny$bfb_event = phylogeny$distance > eps
  phylogeny$birth_time = NA
  phylogeny$hotspot_gained = NA
  phylogeny$is_alive = phylogeny$type == "observed"
  phylogeny$bfb_event[is.na(phylogeny$bfb_event)] = FALSE

  # Generate final cell dataframe
  cell_df = find_bfb_only_cell_data(phylogeny)

  list(cell_df = cell_df, pd_list = pd_list)
}
