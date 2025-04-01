# Helper functions for MFNJ ####

custom_distance <- function(v1, v2, alpha, cn_weight) {
  diff = abs(v1 - v2)
  if (any(diff%%2 == 0 & diff > 0)) {
    diff[diff%%2 == 0 & diff > 0] = alpha
  }
  dist = sqrt(sum((diff)^2))
  #dist = dist - cn_weight * max(c(v1, v2))
  #sum((diff)**2) - sum(diff%%2 == 0 & diff > 0) * cn_weight
  sum((diff)**2) - cn_weight * max(c(v1, v2))
  #dist
}

# custom_distance <- function(v1, v2, alpha, cn_weight) {
#   d1 = sqrt(sum((v1 - v2)^2))
#   normal = rep(1, length(v1))
#   vmerged = (v2 + v2) / 2
#   d2 = sqrt(sum((vmerged - normal)^2))
#   d1 * alpha + (1 - alpha) * d2
# }

euclidean_distance = function(v1, v2) {
  sqrt(sum((v1 - v2)^2))
}

# Function to test detect BFB compatible cells
bfb_compatibility = function(v1, v2, cn_weight) {
  diff = abs(v1 - v2)
  bfb_compatible_differences = diff[diff%%2 == 0 & diff > 0]
  len_bfb_compatible_differences = length(bfb_compatible_differences)
  if (len_bfb_compatible_differences == 0) {
    1
  } else {
    bfb_compatible_percentage = len_bfb_compatible_differences / length(v1)
    1 - bfb_compatible_percentage
  }
}

compute_distance_matrix = function(vectors, alpha = 1/2, cn_weight = 0) {
  n = nrow(vectors)
  if (n == 1) {
    return(matrix(0, nrow = n, ncol = n))
  }

  dist_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist_matrix[j, i] = dist_matrix[i, j] = custom_distance(vectors[i, ], vectors[j, ], alpha, cn_weight)
    }
  }
  dist_matrix
}

compute_compatibility_matrix = function(vectors, cn_weight) {
  n = nrow(vectors)
  if (n == 1) {
    return(matrix(0, nrow = n, ncol = n))
  }

  m <- matrix(1, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      m[j, i] = m[i, j] = bfb_compatibility(vectors[i, ], vectors[j, ], cn_weight = cn_weight)
    }
  }
  m
}

# Distance functions (optimized versions)
node_to_node_distance <- function(i, j, D) D[i, j]

node_to_OTU_distance <- function(node, OTU, D) {
  distances <- D[node, OTU]
  intra_distances <- D[OTU, OTU][upper.tri(D[OTU, OTU])]
  n <- length(OTU)
  (sum(distances)/n) - (sum(intra_distances)/(n*(n-1)))
}

OTU_to_OTU_distance <- function(u, v, D) {
  inter_dist <- sum(D[u, v])
  intra_u <- sum(D[u, u][upper.tri(D[u, u])])
  intra_v <- sum(D[v, v][upper.tri(D[v, v])])
  n_u <- length(u)
  n_v <- length(v)

  (inter_dist/(n_u*n_v)) - (intra_u/(n_u*(n_u-1))) - (intra_v/(n_v*(n_v-1)))
}

# Iterative merging function
# iterative_merge <- function(D, min_size = 2) {
#   otus_history = list()
#   iteration = 1
#   while(nrow(D) > min_size) {
#     new_OTUs <- find_new_OTUs(D)
#
#     # If no merges found, break to prevent infinite loop
#     if(length(new_OTUs) == 0) break
#
#     # Build new distance matrix
#     all_nodes <- colnames(D)
#     single_nodes <- setdiff(all_nodes, unlist(new_OTUs))
#     elements <- c(new_OTUs, as.list(single_nodes))
#
#     # Create new matrix
#     new_D <- matrix(0, nrow = length(elements), ncol = length(elements))
#     rownames(new_D) <- colnames(new_D) <- sapply(elements, function(x) paste(x, collapse = "_"))
#
#     # Compute distances
#     for (i in 1:(nrow(new_D)-1)) {
#       for (j in (i+1):ncol(new_D)) {
#         x <- elements[[i]]
#         y <- elements[[j]]
#
#         if (length(x) == 1 && length(y) == 1) {
#           dist <- node_to_node_distance(x, y, D)
#         } else if (length(x) == 1) {
#           dist <- node_to_OTU_distance(x, y, D)
#         } else if (length(y) == 1) {
#           dist <- node_to_OTU_distance(y, x, D)
#         } else {
#           dist <- OTU_to_OTU_distance(x, y, D)
#         }
#
#         new_D[i, j] <- new_D[j, i] <- dist
#       }
#     }
#
#     otus_history[[iteration]] = new_OTUs
#     D <- new_D
#     cat("Merged to", nrow(D), "clusters\n")
#     iteration = iteration + 1
#   }
#   return(list(D=D, otus_history=otus_history))
# }

build_new_pseudocell = function(new_otu, pseudo_cells) {
  pseudo_cell = matrix(colMeans(pseudo_cells[new_otu,]), nrow = 1)

  #lapply(new_otu, function(u) {gsub("\\|", "", u)}) %>% unlist()
  new_name = paste0(lapply(new_otu, function(u) {gsub("\\|", "", u)}) %>% unlist(), collapse = "-")
  new_name = paste0(unlist(strsplit(new_name, "-")) %>% sort(), collapse = "-")
  new_name = paste0("|",new_name,"|")
  #new_name = paste0("|", new_name, "|")
  rownames(pseudo_cell) = new_name
  pseudo_cell
}

filter_pseudo_cells = function(pseudo_cells, nu) {
  good_idx = which(!rownames(pseudo_cells) %in% unlist(nu))

  filtered_pseudo_cell = matrix(pseudo_cells[good_idx,], nrow = length(good_idx), ncol = ncol(pseudo_cells))
  rownames(filtered_pseudo_cell) = rownames(pseudo_cells)[good_idx]
  filtered_pseudo_cell
}

get_tree_metrics <- function(tree1_text, tree2_text) {
  tree1 <- ape::read.tree(text = tree1_text)
  tree2 <- ape::read.tree(text = tree2_text)

  common_taxa <- intersect(tree1$tip.label, tree2$tip.label)
  tree1 <- ape::drop.tip(tree1, setdiff(tree1$tip.label, common_taxa))
  tree2 <- ape::drop.tip(tree2, setdiff(tree2$tip.label, common_taxa))

  rf_dist <- phangorn::RF.dist(tree1, tree2)
  #cat("Robinson-Foulds Distance:", rf_dist, "\n")

  sq_status <- Quartet::QuartetStatus(tree1, tree2)
  sim_metrics <- Quartet::SimilarityMetrics(sq_status)
  #cat("Sim metrics:", sim_metrics, "\n")

  # spr_dist <- SPR.dist(tree1, tree2)
  # cat("SPR Distance:", spr_dist, "\n")

  names = c("Robinson-Foulds Distance", colnames(sim_metrics))
  values = c(rf_dist, as.numeric(unlist(sim_metrics)))

  dplyr::tibble(metric = names, value = values)
}

tibble_to_newick <- function(cell_data) {
  # Helper function to recursively build the tree
  build_tree <- function(node) {
    # Find children of the current node
    children <- cell_data %>% dplyr::filter(.data$parent_id == node) %>% dplyr::pull(.data$cell_id)

    if (length(children) == 0) {
      # If no children, return the node itself
      return(node)
    } else {
      # Recursively build subtrees for each child
      subtree <- paste(sapply(children, build_tree), collapse = ",")
      return(paste0("(", subtree, ")", node))
    }
  }

  # Identify the root node (cells with no parent)
  root <- cell_data %>% dplyr::filter(is.na(.data$parent_id)) %>% dplyr::pull(.data$cell_id)

  if (length(root) != 1) {
    stop("Error: There must be exactly one root node.")
  }

  # Build the tree starting from the root
  newick_tree <- paste0(build_tree(root), ";")

  return(newick_tree)
}

# process_mfnj_results = function(mfnj_res, eps=.01)  {
#   pseudo_cells_history = mfnj_res$pseudo_cells_history
#   pseudo_cells_df_history = mfnj_res$pseudo_cells_df_history
#   # Retrieve cell_data from psuedo_cells_history
#   all_cells = pseudo_cells_df_history[[1]] %>% dplyr::pull(cell_id)
#   #all_cells = paste0(unlist(pseudo_cells_history), collapse = "|")
#   #all_cells = unique(unlist(strsplit(x = all_cells, split = "-", fixed = T)))
#   cell_df = dplyr::tibble(cell_id = all_cells, parent_id = NA, type = "observed")
#
#   for (i in 1:length(pseudo_cells_history)) {
#     pseudos = pseudo_cells_history[[i]]
#     for (ps in pseudos) {
#       cell_df = cell_df %>%
#         dplyr::rowwise() %>%
#         dplyr::mutate(parent_id = ifelse(is.na(parent_id) & grepl(gsub("\\|", "-", cell_id), gsub("\\|", "-", ps), fixed = TRUE), ps, parent_id)) %>%
#         dplyr::ungroup()
#       cell_df = dplyr::bind_rows(
#         cell_df,
#         dplyr::tibble(cell_id = ps, parent_id = NA, type = "pseudo")
#       )
#     }
#   }
#
#   # Compute distance between every cell and their parent
#   unique_pseudo_cells = pseudo_cells_df_history %>% do.call("bind_rows", .) %>% dplyr::select(!iter) %>% dplyr::distinct()
#   cell_df$distance = lapply(1:nrow(cell_df), function(i) {
#     cell_and_parent_ids = c(cell_df[i,]$cell_id, cell_df[i,]$parent_id)
#     cells_mat = unique_pseudo_cells %>%
#       dplyr::filter(cell_id %in% cell_and_parent_ids) %>%
#       dplyr::select(!cell_id) %>%
#       as.matrix()
#     if (nrow(cells_mat) == 2) return(euclidean_distance(cells_mat[1,], cells_mat[2,]))
#     return(NA)
#   }) %>% unlist()
#
#   cell_df$distance_to_normal = lapply(1:nrow(cell_df), function(i) {
#     target = cell_df[i,]$cell_id
#     cells_mat = unique_pseudo_cells %>%
#       dplyr::filter(cell_id == target) %>%
#       dplyr::select(!cell_id) %>%
#       as.matrix()
#     if (nrow(cells_mat) == 1) return(euclidean_distance(cells_mat[1,], rep(1, ncol(cells_mat))))
#     return(NA)
#   }) %>% unlist()
#
#   cell_df = cell_df %>% dplyr::mutate(bfb_event = distance > eps)
#
#   unique_pseudo_cells = c(cell_df$cell_id, cell_df$parent_id)[grepl("|", c(cell_df$cell_id, cell_df$parent_id), fixed = TRUE)]
#   pseudo_cells_id_mapping = structure(
#     as.numeric(factor(unique_pseudo_cells)),
#     names = unique_pseudo_cells
#   )
#
#   cell_df$cell_id = lapply(1:nrow(cell_df), function(i) {
#     if (cell_df[i,]$type == "pseudo") {
#       paste0("p", pseudo_cells_id_mapping[[cell_df[i,]$cell_id]])
#     } else {
#       cell_df[i,]$cell_id
#     }
#   }) %>% unlist()
#
#   cell_df$parent_id = lapply(1:nrow(cell_df), function(i) {
#     if (grepl("|", cell_df[i,]$parent_id, fixed = TRUE)) {
#       paste0("p", pseudo_cells_id_mapping[[cell_df[i,]$parent_id]])
#     } else {
#       cell_df[i,]$parent_id
#     }
#   }) %>% unlist()
#
#   cell_df$birth_time = NA
#   cell_df$hotspot_gained =  NA
#   cell_df$is_alive = cell_df$type == "observed"
#
#   cell_df
# }

# Function to extract unique cell references
divide_cell_groups <- function(cell_string) {
  # Remove leading and trailing ||
  cell_string <- gsub("^\\|\\||\\|\\|$", "", cell_string)

  # Split the string at repetitions of ||
  # This will preserve subgroups
  groups <- strsplit(cell_string, "\\|\\|\\|\\|")[[1]]

  # Wrap each group with vertical bars, ensuring || are preserved within subgroups
  groups <- paste0("||", groups, "||")

  return(groups)
}

# Function to compute containment matrix
compute_containment_matrix <- function(input_matrix) {
  # Extract unique cells
  input_cells <- rownames(input_matrix)
  groups = lapply(input_cells, divide_cell_groups)

  # Create an empty matrix
  n <- length(input_cells)
  containment_matrix <- matrix(FALSE, nrow = n, ncol = n)
  rownames(containment_matrix) <- colnames(containment_matrix) <- input_cells
  diag(containment_matrix) = TRUE

  i = 31
  j = 27

  # Populate the matrix
  for (i in 1:n) {
    for (j in 1:n) {
      name_i = input_cells[i]
      name_j = input_cells[j]
      # Check if cells[j] is contained in input_string
      first_check = grepl(name_j, name_i, fixed = TRUE) | grepl(name_i, name_j, fixed = TRUE)
      second_check = lapply(groups, function(group) {
        any(grepl(name_i, group, fixed = TRUE)) & any(grepl(name_j, group, fixed = TRUE))
      }) %>% unlist() %>% any()
      containment_matrix[i, j] = first_check | second_check
    }
  }

  containment_matrix
}

# create_pairwise_mask <- function(input_matrix) {
#   # Get the rownames
#   rownames_list <- rownames(input_matrix)
#   state_sets <- lapply(rownames_list, function(s) unique(unlist(strsplit(s, " "))))
#
#   # Compute the pairwise containment matrix (symmetric)
#   n <- length(state_sets)
#   containment_matrix <- matrix(FALSE, n, n)
#
#   for (i in 1:n) {
#     for (j in i:n) {  # Start from i to avoid redundant checks
#       if (all(state_sets[[i]] %in% state_sets[[j]]) || all(state_sets[[j]] %in% state_sets[[i]])) {
#         containment_matrix[i, j] <- TRUE
#         containment_matrix[j, i] <- TRUE  # Ensure symmetry
#       }
#     }
#   }
# }


# find_new_OTU_using_compatiblity_score_and_masking = function(pseudo_cells, alpha, cn_weight) {
#   mask = compute_containment_matrix(pseudo_cells)
#   D = compute_distance_matrix(pseudo_cells, alpha = alpha, cn_weight = cn_weight)
#   #D_eps = compute_distance_matrix(pseudo_cells, alpha = 1e-4, cn_weight = 0)
#   D_eps = compute_distance_matrix(pseudo_cells, alpha = 1e-4, cn_weight = 1e-4)
#   diag(D_eps) = diag(D) = Inf
#   rownames(D_eps) = colnames(D_eps) = rownames(D) = colnames(D) = rownames(mask)
#   D[mask] = Inf
#   D_eps[mask] = Inf
#
#   cell_to_print = c("|cell_45|", "|cell_46|")
#   print(D[rownames(D) %in% cell_to_print, rownames(D) %in% cell_to_print])
#   print(D_eps[rownames(D_eps) %in% cell_to_print, rownames(D_eps) %in% cell_to_print])
#
#   temp_D = D + D_eps
#   min_dist = min(temp_D)
#
#   print(paste0("min_distance = ", min_dist))
#
#   pairs <- which(temp_D == min_dist & lower.tri(temp_D), arr.ind = TRUE)
#   pairs <- cbind(rownames(temp_D)[pairs[,1]], colnames(temp_D)[pairs[,2]])
#
#   # Convert to list of named sets
#   pair_sets <- apply(pairs, 1, function(x) list(sort(x)))
#   pair_sets <- lapply(pair_sets, function(x) x[[1]])
#
#   # Merge overlapping sets
#   if (length(pair_sets) > 0) {
#     changed <- TRUE
#     while (changed) {
#       changed <- FALSE
#       new_sets <- list()
#
#       for (current in pair_sets) {
#         merged <- FALSE
#         for (i in seq_along(new_sets)) {
#           if (length(intersect(current, new_sets[[i]]))) {
#             new_sets[[i]] <- union(new_sets[[i]], current)
#             merged <- changed <- TRUE
#             break
#           }
#         }
#         if (!merged) new_sets <- c(new_sets, list(current))
#       }
#       pair_sets <- new_sets
#     }
#   }
#
#   # Obtain only maximal vector in each pair set
#   lapply(pair_sets, extract_maximal_cells)
# }


extract_maximal_cells <- function(input_strings) {
  # Extract all individual cell references (e.g., "cell_25", "cell_42", etc.)
  all_cells <- unique(unlist(regmatches(input_strings, gregexpr("cell_\\d+", input_strings))))

  # Initialize a list to store maximal cells
  maximal_cells <- character(0)

  for (s in input_strings) {
    # Check if this string contains any other cell references besides its own
    cells_in_s <- unlist(regmatches(s, gregexpr("cell_\\d+", s)))
    other_cells <- setdiff(all_cells, cells_in_s)

    # If this string does not fully contain any other cell outside its own, keep it
    is_maximal <- TRUE
    for (other_s in input_strings) {
      if (s == other_s) next  # Skip self-comparison
      if (grepl(gsub("\\|", "\\\\|", s), other_s, fixed = FALSE)) {
        is_maximal <- FALSE
        break
      }
    }

    if (is_maximal) {
      maximal_cells <- c(maximal_cells, s)
    }
  }

  # Remove duplicates (if any)
  maximal_cells <- unique(maximal_cells)

  return(maximal_cells)
}

